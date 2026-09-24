#' Generate a MultiMAP analysis report template (binomial endpoint)
#'
#' Writes a Quarto (.qmd) or R Markdown (.Rmd) report template for a binomial
#' MultiMAP object (the output of \code{\link{SCI_binom}} or
#' \code{\link{PI_binom}}), suitable as a starting point for the statistical
#' analysis plan (SAP) of a clinical study. The fitted object is saved next to
#' the report as an \code{.rds} file; the template loads it and rebuilds every
#' table and figure with live code chunks, so the document re-renders
#' reproducibly. The italicised prose is placeholder guidance for the user to
#' replace.
#'
#' Both object types share a block of three sections: the prior (\code{$prior})
#' as a table and density plot; the prediction interval (\code{$pred_int}) for
#' the concurrent control proportion as a table; and the robustified prior
#' (\code{$prior_rob}) as a table and density plot. \code{\link{SCI_binom}}
#' output additionally gets an analysis overview, the historical control data,
#' the current trial, the posterior distributions, the simultaneous credible
#' intervals (as a table, a forest plot and an interpretation), and a reference
#' list. \code{\link{PI_binom}} output gets the overview, the historical control
#' data, the shared block and the reference list.
#'
#' The prior / robustified-prior / posterior mixture densities are drawn with
#' RBesT's own \code{plot} method; the prediction-interval figure and the
#' simultaneous-credible-interval figure use the generic
#' \code{\link{plot.MultiMAP}} method (\code{which = 1} and \code{which = 2}).
#'
#' For a multi-arm borrowing object (\code{histdat} supplied to
#' \code{\link{SCI_binom}} / \code{\link{PI_binom}} as a named list), the template
#' instead contains one \dQuote{Borrowed arms} section per borrowed arm -- each
#' with that arm's historical data, prior, prediction interval, robustified prior
#' and (for \code{SCI_binom} input) posterior -- followed by the non-borrowed
#' groups and the simultaneous credible intervals.
#'
#' @param x a binomial MultiMAP object produced by \code{\link{SCI_binom}} or
#'   \code{\link{PI_binom}} (single-arm or multi-arm).
#' @param path output directory. \code{NULL} (default) writes into
#'   \code{tempdir()} (CRAN policy: do not write to the user's filespace unless a
#'   location is given explicitly); pass an explicit directory (e.g.
#'   \code{getwd()}) to write elsewhere.
#' @param file file name of the report within \code{path}. \code{NULL} (default)
#'   uses \code{"MultiMAP_report.qmd"} (qmd) or \code{"MultiMAP_report.Rmd"}
#'   (rmd). A full path may also be given, in which case its directory overrides
#'   \code{path}.
#' @param title report title (YAML).
#' @param author report author (YAML).
#' @param format \code{"qmd"} (Quarto, default) or \code{"rmd"} (R Markdown).
#' @param overwrite overwrite an existing report file? Default \code{FALSE}.
#'
#' @return The report file path, invisibly. As a side effect it writes the report
#'   file and the \code{"<file>_object.rds"} it loads.
#'
#' @seealso \code{\link{SCI_binom}}, \code{\link{PI_binom}}, \code{\link{report_norm}}
#'
#' @examples
#' \donttest{
#' histdat <- data.frame(events     = c(11, 8, 14, 9),
#'                       non_events = c(29, 32, 26, 31))
#' newdat  <- data.frame(group      = c("control", "d1", "d454"),
#'                       events     = c(8, 14, 17),
#'                       non_events = c(22, 16, 13))
#' res <- SCI_binom(histdat, newdat, prior_est = "MAP",
#'                  contr = "pi_diff", seed = 1)
#' report_binom(res)   # writes into tempdir()
#' }
#' @export
report_binom <- function(x,
                    path      = NULL,
                    file      = NULL,
                    title     = "MultiMAP analysis report",
                    author    = "",
                    format    = c("qmd", "rmd"),
                    overwrite = FALSE){

        # Validate: binomial MultiMAP object from SCI_binom() / PI_binom()
        if(!inherits(x, "MultiMAP")){
                stop("x is not of class MultiMAP")
        }
        if(!inherits(x, "MultiMAP_binom")){
                stop("report_binom() expects a binomial MultiMAP object (from SCI_binom() or PI_binom()); use report_norm() for normal endpoints.")
        }
        # SCI_binom() objects carry class "SCI"; PI_binom() objects carry class
        # "PI" and provide only the prior / prediction interval.
        is_sci <- inherits(x, "SCI")

        # Multi-arm borrowing objects carry $borrow_arms and store their prior /
        # posterior / prediction-interval components as per-arm named lists.
        multi  <- !is.null(x$borrow_arms)

        format <- match.arg(format)

        # Output directory: default to a session temporary directory rather than
        # the working directory (CRAN policy: do not write to the user's
        # filespace unless a location is given explicitly).
        if(is.null(path)){
                path <- tempdir()
        }
        if(!is.character(path) || length(path) != 1){
                stop("'path' must be a single directory path or NULL.")
        }

        # File name within 'path'. A NULL 'file' uses the default name; a bare
        # file name is placed inside 'path'; a file given with its own directory
        # keeps that directory (overriding 'path').
        if(is.null(file)){
                file <- if(format == "qmd") "MultiMAP_report.qmd" else "MultiMAP_report.Rmd"
        }
        if(identical(dirname(file), ".")){
                file <- file.path(path, file)
        }
        if(file.exists(file) && !isTRUE(overwrite)){
                stop("'", file, "' already exists; use overwrite = TRUE to replace it.")
        }

        # Save the fitted object beside the report; the template loads it by name.
        data_file <- paste0(tools::file_path_sans_ext(file), "_object.rds")
        saveRDS(x, data_file)

        #-----------------------------------------------------------------------
        # YAML front matter (format-specific)

        if(format == "qmd"){
                yaml <- c('---',
                          sprintf('title: "%s"', title),
                          sprintf('author: "%s"', author),
                          'date: today',
                          'format:',
                          '  html:',
                          '    toc: true',
                          '    number-sections: true',
                          '---', '')
        } else {
                yaml <- c('---',
                          sprintf('title: "%s"', title),
                          sprintf('author: "%s"', author),
                          'date: "`r Sys.Date()`"',
                          'output:',
                          '  html_document:',
                          '    toc: true',
                          '    number_sections: true',
                          '---', '')
        }

        #-----------------------------------------------------------------------
        # Body (identical for qmd / rmd). A raw string keeps the embedded R code
        # (backslashes, backticks, chunk options) verbatim; {{DATA_FILE}} is the
        # only token substituted below.

        preamble <- r"---(*This is an auto-generated MultiMAP report template. Replace the italicised
guidance with study-specific text; the tables and figures are produced from the
fitted object and update whenever the document is re-rendered.*

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE, message = FALSE, warning = FALSE,
                      fig.width = 7, fig.height = 4)
library(RBesT)
library(ggplot2)
library(knitr)

# MultiMAP provides SCI_binom() / PI_binom() and the summary / get_output
# methods used below.
library(MultiMAP)

# Fitted object produced by SCI_binom() / PI_binom() and saved by report_binom():
x <- readRDS("{{DATA_FILE}}")

# Render an RBesT beta mixture (rows w, a, b) as a tidy table, with the implied
# component mean a / (a + b).
mix_table <- function(m){
  cn <- colnames(m); if(is.null(cn)) cn <- paste0("comp", seq_len(ncol(m)))
  data.frame(component = cn,
             weight = round(m["w", ], 3),
             a      = round(m["a", ], 3),
             b      = round(m["b", ], 3),
             mean   = round(m["a", ] / (m["a", ] + m["b", ]), 3),
             row.names = NULL)
}

# Plot a mixture as a density. automixfit() tags the MAP prior with the extra
# class "EM", which would dispatch to a fit diagnostic; strip it so plot() shows
# the mixture density (RBesT::plot.mix).
as_mix <- function(m){
  keep <- intersect(class(m), c("normMix", "betaMix", "gammaMix", "mix"))
  class(m) <- keep
  m
}

# Coverage level and null value for the contrasts (ratios test against 1,
# differences against 0).
conf     <- if (is.null(x$alpha)) NA else 1 - x$alpha
null_val <- if (identical(x$contr, "pi_diff")) 0 else 1

# Weight of the non-informative (robust) mixture component: prior vs posterior
w_ni_prior <- if (is.data.frame(x$weights)) x$weights$w_prior[2] else NA
w_ni_post  <- if (is.data.frame(x$weights)) x$weights$w_post[2]  else NA

# Weight of the robust (non-informative) component, recoverable for any object
# as the widest (least-concentrated, smallest a + b) component of the
# robustified prior.
w_rob <- unname(x$prior_rob["w", which.min(x$prior_rob["a", ] + x$prior_rob["b", ])])
```
)---"

        #--- type-specific analysis overview -----------------------------------
        overview_sci <- r"---(
# Analysis overview

*Describe the endpoint, the borrowing objective and the trial design. This
analysis borrows historical control information into the concurrent control
group (CCG) of the current trial via a robustified prior for the control
response proportion (Schmidli et al., 2014), and compares each treatment group
with the concurrent control using simultaneous credible intervals at the
`r 100 * conf`% level. The prior is derived with the RBesT package (Weber et al.,
2021). The simultaneous credible intervals are computed from a joint sample of
the posterior distributions using the rank-based method of Besag et al. (1995).*
)---"

        overview_pi <- r"---(
# Analysis overview

*Describe the endpoint and the borrowing objective. This analysis summarises the
historical control information as a robustified prior for the control response
proportion (Schmidli et al., 2014), derived with the RBesT package (Weber et
al., 2021), and reports the resulting `r 100 * conf`% prediction interval for the
proportion of a future concurrent control group.*
)---"

        #--- historical control data (shared) ----------------------------------
        historical <- r"---(
# Historical control data

*List the historical studies that contribute control information to the prior.*

```{r}
kable(x$histdat,
      caption = "Historical control studies (events, non-events).")
```
)---"

        #--- current trial (SCI_binom() only) ----------------------------------
        current_trial <- r"---(
# Current trial

*Describe the arms of the current trial. Row 1 is the concurrent control (CCG);
the remaining rows are the treatment groups.*

```{r}
nd   <- x$newdat
show <- nd[, 1:3]
show$proportion <- round(nd$mean, 3)
kable(show,
      caption = "Current trial: concurrent control (row 1) and treatment groups (events, non-events, observed proportion).")
```
)---"

        #--- shared block: prior, prediction interval, robustified prior --------
        shared_map <- r"---(
# Prior

*The prior summarises the historical control information as a beta mixture; it
is the informative prior for the concurrent control proportion.*

```{r}
kable(mix_table(x$prior),
      caption = "Prior for the concurrent control proportion (beta mixture components).")
```

```{r}
#| fig-cap: "Prior density for the concurrent control proportion (black: mixture; dashed: components)."
plot(as_mix(x$prior))
```

# Prediction interval

*The `r 100 * conf`% prediction interval for the proportion of a concurrent
control group, derived from the prior-predictive (beta-binomial) distribution.*

```{r}
pint <- x$pred_int
kable(data.frame(median = round(unname(pint["median"]), 3),
                 lower  = round(unname(pint["lower"]),  3),
                 upper  = round(unname(pint["upper"]),  3)),
      caption = sprintf("%.0f%% prediction interval for the concurrent control proportion.", 100 * conf))
```

```{r}
#| fig-cap: "Beta-binomial predictive on the proportion scale; bars inside the prediction interval are highlighted. For SCI_binom() input the observed concurrent control proportion is drawn as a vertical line -- black if covered, red if not."
plot(x, which = 1)
```

# Robustified prior

*The prior is robustified by mixing in a weakly-informative (non-informative)
component with a weight of `r round(100 * w_rob, 1)`%. This guards against
prior-data conflict: if the concurrent control disagrees with the historical
data, the informative component is automatically down-weighted.*

```{r}
kable(mix_table(x$prior_rob),
      caption = "Robustified prior for the concurrent control proportion (beta mixture components).")
```

```{r}
#| fig-cap: "Robustified prior density (black: mixture; dashed: components)."
plot(as_mix(x$prior_rob))
```
)---"

        #--- posterior + simultaneous credible intervals (SCI_binom() only) ----
        posterior_sci <- r"---(
# Posterior distributions

## Concurrent control

*Posterior of the concurrent control proportion, obtained by updating the
robustified prior with the current control data. In the posterior, the
non-informative component has a weight of `r round(100 * w_ni_post, 1)`% (prior weight `r round(100 * w_ni_prior, 1)`%).
It is thus `r if (isTRUE(w_ni_post > w_ni_prior)) "up-weighted" else "down-weighted"` relative to its prior weight,
`r if (isTRUE(w_ni_post > w_ni_prior)) "which indicates prior-data conflict and hence reduced borrowing of historical information" else "which indicates agreement between the concurrent control and the historical data, so the borrowed information is largely retained"`.*

```{r}
kable(mix_table(x$posterior),
      caption = "Posterior of the concurrent control proportion (beta mixture components).")
```

```{r}
#| fig-cap: "Posterior density of the concurrent control proportion."
plot(as_mix(x$posterior))
```

## Treatment groups

*Each treatment posterior is Beta(1 + events, 1 + non-events), i.e. a uniform
Beta(1, 1) prior updated with that arm's data; no historical information is
borrowed for the treatment arms.*

```{r}
trt <- x$newdat[-1, , drop = FALSE]
kable(data.frame(group     = trt[[1]],
                 post_a    = trt$posta,
                 post_b    = trt$postb,
                 post_mean = round(trt$posta / (trt$posta + trt$postb), 3)),
      caption = "Treatment-group posteriors (beta).")
```

# Simultaneous credible intervals

*Treatment-vs-control contrasts with simultaneous credible intervals at the
`r 100 * conf`% level.*

```{r}
kable(x$comp,
      caption = "Simultaneous credible intervals for the treatment-vs-control contrasts.")
```

```{r}
#| fig-cap: "Simultaneous credible intervals; the dashed line marks the null value."
plot(x, which = 2)
```

## Interpretation

```{r}
#| results: asis
cmp    <- x$comp
lab    <- cmp[[1]]
covers <- cmp$lower.ci <= null_val & cmp$upper.ci >= null_val
contr_word <- switch(as.character(x$contr),
                     RR = "ratio of proportions",
                     OR      = "odds ratio",
                     pi_diff = "difference of proportions",
                     "contrast")
cat(sprintf("The null value for a %s is **%g**.\n\n", contr_word, null_val))
for (i in seq_len(nrow(cmp))) {
  cat(sprintf("- **%s**: [%.3g, %.3g] %s the null value -- %s.\n",
              lab[i], cmp$lower.ci[i], cmp$upper.ci[i],
              if (covers[i]) "covers" else "excludes",
              if (covers[i]) "no credible treatment effect"
                        else "a credible treatment effect"))
}
```

*A contrast whose simultaneous credible interval excludes the null value provides
credible evidence of a treatment-vs-control effect at the chosen simultaneous
level; a contrast covering the null does not.*
)---"

        #--- references ---------------------------------------------------------
        refs_sci <- r"---(
# References

Menssen, M., Kneuer, C., Akyianu, G., Roever, C., Friede, T., & Schaarschmidt, F.
(2026). Including historical control data in simultaneous inference for pre-clinical
multi-arm studies. *arXiv preprint* arXiv:2603.11730.
https://doi.org/10.48550/arXiv.2603.11730

Besag, J., Green, P., Higdon, D., & Mengersen, K. (1995). Bayesian computation and
stochastic systems. *Statistical Science, 10*(1), 3--66.
http://www.jstor.org/stable/2246224

Schmidli, H., Gsteiger, S., Roychoudhury, S., O'Hagan, A., Spiegelhalter, D., &
Neuenschwander, B. (2014). Robust meta-analytic-predictive priors in clinical trials
with historical control information. *Biometrics, 70*(4), 1023--1032.
https://doi.org/10.1111/biom.12242

Weber, S., Li, Y., Seaman, J. W., Kakizume, T., & Schmidli, H. (2021). Applying
meta-analytic-predictive priors with the R Bayesian evidence synthesis tools.
*Journal of Statistical Software, 100*(19), 1--32.
https://doi.org/10.18637/jss.v100.i19
)---"

        refs_pi <- r"---(
# References

Menssen, M., Kneuer, C., Akyianu, G., Roever, C., Friede, T., & Schaarschmidt, F.
(2026). Including historical control data in simultaneous inference for pre-clinical
multi-arm studies. *arXiv preprint* arXiv:2603.11730.
https://doi.org/10.48550/arXiv.2603.11730

Schmidli, H., Gsteiger, S., Roychoudhury, S., O'Hagan, A., Spiegelhalter, D., &
Neuenschwander, B. (2014). Robust meta-analytic-predictive priors in clinical trials
with historical control information. *Biometrics, 70*(4), 1023--1032.
https://doi.org/10.1111/biom.12242

Weber, S., Li, Y., Seaman, J. W., Kakizume, T., & Schmidli, H. (2021). Applying
meta-analytic-predictive priors with the R Bayesian evidence synthesis tools.
*Journal of Statistical Software, 100*(19), 1--32.
https://doi.org/10.18637/jss.v100.i19
)---"

        #=======================================================================
        # Multi-arm variant: one section per borrowed arm (each with its own
        # historical data, prior, prediction interval, robustified prior and --
        # for SCI_binom() input -- posterior).
        #=======================================================================

        preamble_multi <- r"---(*This is an auto-generated MultiMAP multi-arm report template. Replace the
italicised guidance with study-specific text; the tables and figures are
produced from the fitted object and update whenever the document is re-rendered.*

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE, message = FALSE, warning = FALSE,
                      fig.width = 7, fig.height = 4)
library(RBesT)
library(ggplot2)
library(knitr)
library(MultiMAP)

# Fitted multi-arm object produced by SCI_binom() / PI_binom():
x <- readRDS("{{DATA_FILE}}")

# Render an RBesT beta mixture (rows w, a, b) as a tidy table.
mix_table <- function(m){
  cn <- colnames(m); if(is.null(cn)) cn <- paste0("comp", seq_len(ncol(m)))
  data.frame(component = cn,
             weight = round(m["w", ], 3),
             a      = round(m["a", ], 3),
             b      = round(m["b", ], 3),
             mean   = round(m["a", ] / (m["a", ] + m["b", ]), 3),
             row.names = NULL)
}

# Strip the automixfit() "EM" class so plot() shows the mixture density.
as_mix <- function(m){
  keep <- intersect(class(m), c("normMix", "betaMix", "gammaMix", "mix"))
  class(m) <- keep
  m
}

conf     <- if (is.null(x$alpha)) NA else 1 - x$alpha
null_val <- if (identical(x$contr, "pi_diff")) 0 else 1
borrow_arms <- x$borrow_arms
```
)---"

        overview_multi_sci <- r"---(
# Analysis overview

*Describe the endpoint, the borrowing objective and the trial design. This
analysis borrows historical control information into **more than one arm** of the
current trial: each borrowed arm listed below gets its own robustified prior for
the response proportion (Schmidli et al., 2014), updated with that arm's current
data. Treatment groups are compared using simultaneous credible intervals at the
`r 100 * conf`% level (Besag et al., 1995). The priors are derived with the RBesT
package (Weber et al., 2021).*

```{r}
kable(data.frame(arm = borrow_arms), caption = "Arms for which historical information is borrowed.")
```
)---"

        overview_multi_pi <- r"---(
# Analysis overview

*Describe the endpoint and the borrowing objective. This analysis summarises the
historical control information for **more than one arm** as robustified priors
for the response proportion (Schmidli et al., 2014), derived with the RBesT
package (Weber et al., 2021), and reports the resulting `r 100 * conf`% prediction
interval for each borrowed arm.*

```{r}
kable(data.frame(arm = borrow_arms), caption = "Arms for which historical information is borrowed.")
```
)---"

        current_trial_multi <- r"---(
# Current trial

*Describe the arms of the current trial. Borrowed arms use an informative prior;
every other arm uses a non-borrowing Beta(1, 1) prior.*

```{r}
nd   <- x$newdat
show <- nd[, 1:3]
show$proportion <- round(nd$mean, 3)
kable(show, caption = "Current trial groups (events, non-events, observed proportion).")
```
)---"

        arm_block_head <- r"---(
## Arm: <<ARM>>

### Historical control data

```{r}
kable(x$histdat[["<<ARM>>"]],
      caption = "Historical control studies for arm <<ARM>> (events, non-events).")
```

### Prior

```{r}
kable(mix_table(x$prior[["<<ARM>>"]]),
      caption = "Prior for arm <<ARM>> (beta mixture components).")
```

```{r}
#| fig-cap: "Prior density for arm <<ARM>>."
plot(as_mix(x$prior[["<<ARM>>"]]))
```

### Prediction interval

```{r}
pint <- x$pred_int[["<<ARM>>"]]
kable(data.frame(median = round(unname(pint["median"]), 3),
                 lower  = round(unname(pint["lower"]),  3),
                 upper  = round(unname(pint["upper"]),  3)),
      caption = "Prediction interval for arm <<ARM>> (proportion scale).")
```

```{r}
#| fig-cap: "Beta-binomial predictive for arm <<ARM>> on the proportion scale."
plot(x, arm = "<<ARM>>")
```

### Robustified prior

```{r}
kable(mix_table(x$prior_rob[["<<ARM>>"]]),
      caption = "Robustified prior for arm <<ARM>>.")
```

```{r}
#| fig-cap: "Robustified prior density for arm <<ARM>>."
plot(as_mix(x$prior_rob[["<<ARM>>"]]))
```
)---"

        arm_block_post <- r"---(
### Posterior

*Posterior for arm <<ARM>>, obtained by updating its robustified prior with the
arm's current data.*

```{r}
kable(mix_table(x$posterior[["<<ARM>>"]]),
      caption = "Posterior for arm <<ARM>> (beta mixture components).")
```

```{r}
#| fig-cap: "Posterior density for arm <<ARM>>."
plot(as_mix(x$posterior[["<<ARM>>"]]))
```
)---"

        nonborrow_multi <- r"---(
# Non-borrowed groups

*Groups without historical information use a non-borrowing posterior:
Beta(1 + events, 1 + non-events).*

```{r}
allg <- as.character(x$newdat[[1]])
nb   <- setdiff(allg, x$borrow_arms)
if (length(nb)) {
  trt <- x$newdat[match(nb, allg), , drop = FALSE]
  pa  <- 1 + trt[[2]]; pb <- 1 + trt[[3]]
  kable(data.frame(group = trt[[1]], post_a = pa, post_b = pb,
                   post_mean = round(pa / (pa + pb), 3)),
        caption = "Non-borrowed groups: posteriors (beta).")
}
```
)---"

        sci_tail_multi <- r"---(
# Simultaneous credible intervals

*Treatment-vs-control contrasts with simultaneous credible intervals at the
`r 100 * conf`% level.*

```{r}
kable(x$comp,
      caption = "Simultaneous credible intervals for the treatment-vs-control contrasts.")
```

```{r}
#| fig-cap: "Simultaneous credible intervals; the dashed line marks the null value."
plot(x, which = 2)
```

## Interpretation

```{r}
#| results: asis
cmp    <- x$comp
lab    <- cmp[[1]]
covers <- cmp$lower.ci <= null_val & cmp$upper.ci >= null_val
contr_word <- switch(as.character(x$contr),
                     RR = "ratio of proportions",
                     OR      = "odds ratio",
                     pi_diff = "difference of proportions",
                     "contrast")
cat(sprintf("The null value for a %s is **%g**.\n\n", contr_word, null_val))
for (i in seq_len(nrow(cmp))) {
  cat(sprintf("- **%s**: [%.3g, %.3g] %s the null value -- %s.\n",
              lab[i], cmp$lower.ci[i], cmp$upper.ci[i],
              if (covers[i]) "covers" else "excludes",
              if (covers[i]) "no credible treatment effect"
                        else "a credible treatment effect"))
}
```

*A contrast whose simultaneous credible interval excludes the null value provides
credible evidence of a treatment-vs-control effect at the chosen simultaneous
level; a contrast covering the null does not.*
)---"

        #--- assemble -----------------------------------------------------------
        if (multi) {
                arm_tmpl <- if (is_sci) paste0(arm_block_head, arm_block_post) else arm_block_head
                per_arm  <- paste0(
                        vapply(x$borrow_arms,
                               function(a) gsub("<<ARM>>", a, arm_tmpl, fixed = TRUE),
                               character(1)),
                        collapse = "\n")
                per_arm  <- paste0("\n# Borrowed arms\n", per_arm)
                if (is_sci) {
                        body <- paste0(preamble_multi, overview_multi_sci, current_trial_multi,
                                       per_arm, nonborrow_multi, sci_tail_multi, refs_sci)
                } else {
                        body <- paste0(preamble_multi, overview_multi_pi, per_arm, refs_pi)
                }
        } else if (is_sci) {
                body <- paste0(preamble, overview_sci, historical, current_trial,
                               shared_map, posterior_sci, refs_sci)
        } else {
                body <- paste0(preamble, overview_pi, historical, shared_map, refs_pi)
        }
        body <- gsub("{{DATA_FILE}}", basename(data_file), body, fixed = TRUE)

        #-----------------------------------------------------------------------
        # Write the report

        writeLines(c(yaml, body), con = file)

        render_hint <- if(format == "qmd")
                sprintf('quarto::quarto_render("%s")', file)
        else    sprintf('rmarkdown::render("%s")', file)

        message("Report template written to '", file, "' ",
                "(object saved as '", data_file, "').\n",
                "Render it with: ", render_hint)

        invisible(file)
}

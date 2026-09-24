#' Generate a MultiMAP analysis report template
#'
#' Writes a Quarto (.qmd) or R Markdown (.Rmd) report template for a MultiMAP
#' object, suitable as a starting point for the statistical analysis plan (SAP)
#' of a clinical study. The fitted object is saved next to the report as an
#' \code{.rds} file; the template loads it and rebuilds every table and figure
#' with live code chunks, so the document re-renders reproducibly. The
#' italicised prose is placeholder guidance for the user to replace.
#'
#' Both object types share an identical block of three sections (same prose in
#' both reports): the MAP prior (\code{$prior}) as a table and plot; the
#' MAP-based prediction interval (\code{$pred_int}) as a table and plot -- for
#' \code{\link{SCI_norm}} input \code{plot(x)} adds the observed concurrent
#' control mean as a vertical line (black if covered by the interval, red if
#' not); and the robustified MAP prior (\code{$prior_rob}) as a table and plot.
#' \code{SCI_norm} output additionally gets an analysis overview, the historical
#' control data, the current trial, the posterior distributions, the
#' simultaneous credible intervals with an interpretation, and an APA-style
#' reference list. \code{\link{PI_norm}} output gets the overview, the historical
#' control data, the shared block and the reference list.
#'
#' For a multi-arm borrowing object (\code{histdat} supplied to
#' \code{\link{SCI_norm}} / \code{\link{PI_norm}} as a named list), the template
#' instead contains one \dQuote{Borrowed arms} section per borrowed arm -- each
#' with that arm's historical data, MAP prior, prediction interval, robustified
#' prior and (for \code{SCI_norm} input) posterior -- followed by the non-borrowed
#' groups and the simultaneous credible intervals.
#'
#' @param x a MultiMAP object produced by \code{\link{SCI_norm}} or
#'   \code{\link{PI_norm}} (single-arm or multi-arm).
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
#' @seealso \code{\link{SCI_norm}}, \code{\link{PI_norm}}
#'
#' @examples
#' \donttest{
#' histdat <- data.frame(mean = c(9.8, 10.2, 9.5, 10.0),
#'                       sd   = c(2.1, 1.9, 2.3, 2.0),
#'                       n    = c(40, 55, 38, 47))
#' newdat  <- data.frame(group = c("control", "d1", "d454"),
#'                       mean  = c(10.1, 11.4, 12.7),
#'                       sd    = c(2.0, 2.2, 2.1),
#'                       n     = c(30, 30, 30))
#' res <- SCI_norm(histdat, newdat, contr = "mean_diff", seed = 1)
#' report_norm(res)   # writes into tempdir()
#' }
#' @export
report_norm <- function(x,
                   path      = NULL,
                   file      = NULL,
                   title     = "MultiMAP analysis report",
                   author    = "",
                   format    = c("qmd", "rmd"),
                   overwrite = FALSE){

        # Validate: MultiMAP object from SCI_norm() (carries newdat + comp)
        if(!inherits(x, "MultiMAP")){
                stop("x is not of class MultiMAP")
        }
        # SCI_norm() objects carry class "SCI"; PI_norm() objects carry class "PI"
        # and provide only the prior / prediction interval.
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

# MultiMAP provides SCI_norm() / PI_norm() and the plot / print / summary /
# get_output methods used below.
library(MultiMAP)

# Fitted object produced by SCI_norm() and saved by report_norm():
x <- readRDS("{{DATA_FILE}}")

# Render an RBesT normal mixture (rows w, m, s) as a tidy table
mix_table <- function(m){
  cn <- colnames(m); if(is.null(cn)) cn <- paste0("comp", seq_len(ncol(m)))
  data.frame(component = cn,
             weight = round(m["w", ], 3),
             mean   = round(m["m", ], 3),
             sd     = round(m["s", ], 3),
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

# Coverage level and null value for the contrasts
conf     <- if (is.null(x$alpha)) NA else 1 - x$alpha
null_val <- if (identical(x$contr, "mean_ratio")) 1 else 0

# Weight of the non-informative (robust) mixture component: prior vs posterior
w_ni_prior <- if (is.data.frame(x$weights)) x$weights$w_prior[2] else NA
w_ni_post  <- if (is.data.frame(x$weights)) x$weights$w_post[2]  else NA

# Weight of the robust (non-informative) component, recoverable for any object
# as the widest (largest-sd) component of the robustified prior.
w_rob <- unname(x$prior_rob["w", which.max(x$prior_rob["s", ])])
```
)---"

        #--- type-specific analysis overview -----------------------------------
        overview_sci <- r"---(
# Analysis overview

*Describe the endpoint, the borrowing objective and the trial design. This
analysis borrows historical control information into the concurrent control
group (CCG) of the current trial via a robustified meta-analytic-predictive
(MAP) prior (Schmidli et al., 2014), and compares each treatment group with the
concurrent control using simultaneous credible intervals at the `r 100 * conf`%
level. The MAP prior is derived with the RBesT package (Weber et al., 2021).
The simultaneous credible intervals are computed from a joint sample of the
posterior distributions using the rank-based method of Besag et al. (1995).*
)---"

        overview_pi <- r"---(
# Analysis overview

*Describe the endpoint and the borrowing objective. This analysis summarises the
historical control information as a robustified meta-analytic-predictive (MAP)
prior (Schmidli et al., 2014), derived with the RBesT package (Weber et al.,
2021), and reports the resulting `r 100 * conf`% prediction interval for the mean
of a future concurrent control group.*
)---"

        #--- historical control data (shared) ----------------------------------
        historical <- r"---(
# Historical control data

*List the historical studies that contribute control information to the MAP
prior.*

```{r}
kable(x$histdat,
      caption = "Historical control studies (mean, standard deviation, sample size).")
```
)---"

        #--- current trial (SCI_norm() only) -----------------------------------
        current_trial <- r"---(
# Current trial

*Describe the arms of the current trial. Row 1 is the concurrent control (CCG);
the remaining rows are the treatment groups.*

```{r}
cols <- intersect(c("group", "mean", "sd", "n"), colnames(x$newdat))
kable(x$newdat[, cols],
      caption = "Current trial: concurrent control (row 1) and treatment groups.")
```
)---"

        #--- shared block: MAP prior, prediction interval, robustified MAP ------
        # Identical prose in both reports. plot(x) draws the prediction-interval
        # figure: for SCI_norm() input it adds the observed concurrent control
        # mean as a vertical line (black if covered, red if not); for PI_norm()
        # input there is no such line.
        shared_map <- r"---(
# MAP prior

*The MAP prior summarises the historical control information as a normal mixture;
it is the informative prior for the concurrent control mean.*

```{r}
kable(mix_table(x$prior),
      caption = "MAP prior for the concurrent control mean (mixture components).")
```

```{r}
#| fig-cap: "MAP prior density for the concurrent control mean (black: mixture; dashed: components)."
plot(as_mix(x$prior))
```

# Prediction interval

*The `r 100 * conf`% prediction interval for the mean of a concurrent control
group, derived from the MAP prior. When an observed concurrent control mean is
available (SCI_norm() input) it is drawn as a vertical line -- black if it lies
inside the interval, red if it falls outside (indicating prior-data conflict).*

```{r}
pint <- x$pred_int
kable(data.frame(median = round(unname(pint["median"]), 3),
                 lower  = round(unname(pint["lower"]),  3),
                 upper  = round(unname(pint["upper"]),  3)),
      caption = sprintf("%.0f%% prediction interval for the concurrent control mean.", 100 * conf))
```

```{r}
#| fig-cap: "MAP-based predictive density with the prediction interval shaded."
plot(x)
```

# Robustified MAP prior

*The MAP prior is robustified by mixing in a weakly-informative (non-informative)
component with a weight of `r round(100 * w_rob, 1)`%. This guards against
prior-data conflict: if the concurrent control disagrees with the historical
data, the informative component is automatically down-weighted.*

```{r}
kable(mix_table(x$prior_rob),
      caption = "Robustified MAP prior for the concurrent control mean (mixture components).")
```

```{r}
#| fig-cap: "Robustified MAP prior density (black: mixture; dashed: components)."
plot(as_mix(x$prior_rob))
```
)---"

        #--- posterior + simultaneous credible intervals (SCI_norm() only) -----
        posterior_sci <- r"---(
# Posterior distributions

## Concurrent control

*Posterior of the concurrent control mean, obtained by updating the robustified
prior with the current control data. In the posterior, the non-informative
component has a weight of `r round(100 * w_ni_post, 1)`% (prior weight `r round(100 * w_ni_prior, 1)`%).
It is thus `r if (isTRUE(w_ni_post > w_ni_prior)) "up-weighted" else "down-weighted"` relative to its prior weight,
`r if (isTRUE(w_ni_post > w_ni_prior)) "which indicates prior-data conflict and hence reduced borrowing of historical information" else "which indicates agreement between the concurrent control and the historical data, so the borrowed information is largely retained"`.*

```{r}
kable(mix_table(x$posterior),
      caption = "Posterior of the concurrent control mean (mixture components).")
```

```{r}
#| fig-cap: "Posterior density of the concurrent control mean."
plot(as_mix(x$posterior))
```

## Treatment groups

*Each treatment posterior is normal, centred at the observed treatment mean with
standard error sd / sqrt(n) (flat prior on the mean); no historical information
is borrowed for the treatment arms.*

```{r}
trt <- x$newdat[-1, , drop = FALSE]
pm  <- if (!is.null(trt$post_m)) trt$post_m else trt$mean
ps  <- if (!is.null(trt$post_s)) trt$post_s else trt$sd / sqrt(trt$n)
kable(data.frame(group = trt$group, post_mean = round(pm, 3), post_sd = round(ps, 3)),
      caption = "Treatment-group posteriors (normal).")
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
cat(sprintf("The null value for a %s is **%g**.\n\n",
            if (identical(x$contr, "mean_ratio")) "ratio of means" else "difference of means",
            null_val))
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

        #--- references (Besag only for the SCI rank-based method) --------------
        refs_sci <- r"---(
# References

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
        # historical data, MAP prior, prediction interval, robustified prior and
        # -- for SCI_norm() input -- posterior). Built by unrolling a per-arm
        # template so every chunk stays a plain (single-object) chunk.
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

# Fitted multi-arm object produced by SCI_norm() / PI_norm():
x <- readRDS("{{DATA_FILE}}")

# Render an RBesT normal mixture (rows w, m, s) as a tidy table
mix_table <- function(m){
  cn <- colnames(m); if(is.null(cn)) cn <- paste0("comp", seq_len(ncol(m)))
  data.frame(component = cn,
             weight = round(m["w", ], 3),
             mean   = round(m["m", ], 3),
             sd     = round(m["s", ], 3),
             row.names = NULL)
}

# Strip the automixfit() "EM" class so plot() shows the mixture density.
as_mix <- function(m){
  keep <- intersect(class(m), c("normMix", "betaMix", "gammaMix", "mix"))
  class(m) <- keep
  m
}

conf     <- if (is.null(x$alpha)) NA else 1 - x$alpha
null_val <- if (identical(x$contr, "mean_ratio")) 1 else 0
borrow_arms <- x$borrow_arms
```
)---"

        overview_multi_sci <- r"---(
# Analysis overview

*Describe the endpoint, the borrowing objective and the trial design. This
analysis borrows historical control information into **more than one arm** of the
current trial: each borrowed arm listed below gets its own robustified
meta-analytic-predictive (MAP) prior (Schmidli et al., 2014), updated with that
arm's current data. Treatment groups are compared using simultaneous credible
intervals at the `r 100 * conf`% level (Besag et al., 1995). The MAP priors are
derived with the RBesT package (Weber et al., 2021).*

```{r}
kable(data.frame(arm = borrow_arms), caption = "Arms for which historical information is borrowed.")
```
)---"

        overview_multi_pi <- r"---(
# Analysis overview

*Describe the endpoint and the borrowing objective. This analysis summarises the
historical control information for **more than one arm** as robustified
meta-analytic-predictive (MAP) priors (Schmidli et al., 2014), derived with the
RBesT package (Weber et al., 2021), and reports the resulting `r 100 * conf`%
prediction interval for each borrowed arm.*

```{r}
kable(data.frame(arm = borrow_arms), caption = "Arms for which historical information is borrowed.")
```
)---"

        current_trial_multi <- r"---(
# Current trial

*Describe the arms of the current trial. Borrowed arms use a MAP-based prior;
every other arm uses a non-borrowing (flat-prior) posterior.*

```{r}
cols <- intersect(c("group", "mean", "sd", "n"), colnames(x$newdat))
kable(x$newdat[, cols], caption = "Current trial groups.")
```
)---"

        # Per-arm template (SCI: includes the posterior; PI: omitted below).
        arm_block_head <- r"---(
## Arm: <<ARM>>

### Historical control data

```{r}
kable(x$histdat[["<<ARM>>"]],
      caption = "Historical control studies for arm <<ARM>> (mean, sd, n).")
```

### MAP prior

```{r}
kable(mix_table(x$prior[["<<ARM>>"]]),
      caption = "MAP prior for arm <<ARM>> (mixture components).")
```

```{r}
#| fig-cap: "MAP prior density for arm <<ARM>>."
plot(as_mix(x$prior[["<<ARM>>"]]))
```

### Prediction interval

```{r}
pint <- x$pred_int[["<<ARM>>"]]
kable(data.frame(median = round(unname(pint["median"]), 3),
                 lower  = round(unname(pint["lower"]),  3),
                 upper  = round(unname(pint["upper"]),  3)),
      caption = "Prediction interval for arm <<ARM>>.")
```

```{r}
#| fig-cap: "MAP-based predictive density for arm <<ARM>> with the prediction interval shaded."
plot(x, arm = "<<ARM>>")
```

### Robustified MAP prior

```{r}
kable(mix_table(x$prior_rob[["<<ARM>>"]]),
      caption = "Robustified MAP prior for arm <<ARM>>.")
```

```{r}
#| fig-cap: "Robustified MAP prior density for arm <<ARM>>."
plot(as_mix(x$prior_rob[["<<ARM>>"]]))
```
)---"

        arm_block_post <- r"---(
### Posterior

*Posterior for arm <<ARM>>, obtained by updating its robustified prior with the
arm's current data.*

```{r}
kable(mix_table(x$posterior[["<<ARM>>"]]),
      caption = "Posterior for arm <<ARM>> (mixture components).")
```

```{r}
#| fig-cap: "Posterior density for arm <<ARM>>."
plot(as_mix(x$posterior[["<<ARM>>"]]))
```
)---"

        nonborrow_multi <- r"---(
# Non-borrowed groups

*Groups without historical information use a non-borrowing posterior:
normal, centred at the observed mean with standard error sd / sqrt(n).*

```{r}
allg <- as.character(x$newdat$group)
nb   <- setdiff(allg, x$borrow_arms)
if (length(nb)) {
  trt <- x$newdat[match(nb, allg), , drop = FALSE]
  kable(data.frame(group = trt$group,
                   post_mean = round(trt$mean, 3),
                   post_sd   = round(trt$sd / sqrt(trt$n), 3)),
        caption = "Non-borrowed groups: posteriors (normal).")
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
cat(sprintf("The null value for a %s is **%g**.\n\n",
            if (identical(x$contr, "mean_ratio")) "ratio of means" else "difference of means",
            null_val))
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

        #--- assemble: shared blocks are byte-identical in both reports ---------
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

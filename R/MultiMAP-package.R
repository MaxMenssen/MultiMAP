#' MultiMAP: Bayesian dynamic borrowing with simultaneous credible intervals
#'
#' Robustified priors for Bayesian dynamic borrowing from historical control
#' data, with simultaneous credible intervals for treatment-vs-control
#' contrasts. Two endpoint families are supported, each with the same workflow:
#' \itemize{
#'   \item \strong{Normal endpoints} -- prediction intervals
#'     (\code{\link{PI_norm}}), simultaneous credible intervals
#'     (\code{\link{SCI_norm}}), Monte-Carlo operating characteristics
#'     (\code{\link{SIM_norm}}), a data generator (\code{\link{h1_norm}}) and a
#'     report-template generator (\code{\link{report_norm}}). The prior is a
#'     meta-analytic-predictive (MAP) prior via \code{RBesT::gMAP}.
#'   \item \strong{Binomial endpoints} -- the counterparts
#'     \code{\link{PI_binom}}, \code{\link{SCI_binom}}, \code{\link{SIM_binom}},
#'     \code{\link{h1_binom}} and \code{\link{report_binom}}. The prior is either
#'     a MAP prior (\code{RBesT::gMAP}, \code{family = binomial}) or an
#'     empirical-Bayes beta prior estimated from the historical data.
#'   \item \strong{Poisson (count) endpoints} -- the counterparts
#'     \code{\link{PI_pois}}, \code{\link{SCI_pois}}, \code{\link{SIM_pois}},
#'     \code{\link{h1_pois}} and \code{\link{report_pois}} for a gamma-Poisson
#'     rate model. The prior is either a MAP prior (\code{RBesT::gMAP},
#'     \code{family = poisson}) or an empirical-Bayes gamma prior following
#'     Tarone (1982).
#' }
#' Shared \code{plot} / \code{print} / \code{summary} methods and
#' \code{\link{get_output}} work on the objects of both families.
#'
#' @keywords internal
#' @importFrom stats gaussian binomial poisson rnorm rbeta rbinom rgamma rpois plogis qlogis sd median dnorm dgamma sigma
#' @importFrom RBesT gMAP automixfit mixbeta mixgamma robustify postmix rmix ess sigma<-
#' @importFrom utils globalVariables
#' @import ggplot2
"_PACKAGE"

# Silence R CMD check "no visible binding for global variable" NOTEs arising
# from non-standard evaluation:
#   * gMAP() model formula / weights          -> mean, se, study, n
#   * ggplot2 aes() mappings                  -> x, y, density, Distribution,
#                                                comp, estimate, lower.ci, upper.ci,
#                                                prop, pmf, inside
globalVariables(c(
        "mean", "se", "study", "n", "events", "offset",
        "x", "y", "density", "Distribution",
        "comp", "estimate", "lower.ci", "upper.ci",
        "prop", "pmf", "inside", "rate",
        "arm", "width", "xintercept", "col"
))

#-------------------------------------------------------------------------------
# Internal: "informative" x-window for the density / predictive plots.
#
# Returns c(lo, hi): the range where the reference curve `dens` (evaluated at the
# x positions `xv`) is at least `frac` of its peak, always widened to include the
# finite anchor points `include` (e.g. the prediction-interval bounds and the
# observed value) plus a small proportional pad. This makes the x axis focus on
# where the mass / the prediction interval is and drop the uninformative tails
# that heavy-tailed MAP predictives / priors would otherwise stretch across.
.info_window <- function(xv, dens, include = numeric(0), frac = 0.01, pad = 0.05){
        ok <- is.finite(dens) & is.finite(xv)
        if(!any(ok)) return(range(xv[is.finite(xv)]))
        keep <- xv[ok][dens[ok] >= frac * max(dens[ok])]
        lo <- min(keep); hi <- max(keep)
        inc <- include[is.finite(include)]
        if(length(inc)){ lo <- min(lo, inc); hi <- max(hi, inc) }
        w <- hi - lo
        if(!is.finite(w) || w <= 0) w <- max(abs(hi), 1)
        c(lo - pad * w, hi + pad * w)
}

#-------------------------------------------------------------------------------
# Internal: gMAP between-study heterogeneity prior distribution (tau.dist).
#
# The distribution families accepted for the 'tau_dist' argument of the SCI_* /
# PI_* functions, mirroring RBesT::gMAP()'s 'tau.dist'. The built-in default
# 'tau_prior' heuristics are calibrated for the half-normal, so any other family
# must be paired with an explicit 'tau_prior' (as RBesT itself requires).

.tau_dist_choices <- c("HalfNormal", "TruncNormal", "Uniform", "Gamma",
                       "InvGamma", "LogNormal", "TruncCauchy", "Exp", "Fixed")

.check_tau_dist <- function(tau_dist, tau_prior){
        if(!identical(tau_dist, "HalfNormal") && is.null(tau_prior))
                stop("When 'tau_dist' is not \"HalfNormal\", supply a matching 'tau_prior' ",
                     "(the built-in default is only calibrated for the half-normal).")
        invisible(TRUE)
}

# Assemble the $tau_prior output slot: the heterogeneity prior distribution name
# ('dist') plus its parameter matrix ('prior'); NA when no MAP prior was fitted.
.tau_prior_out <- function(tau_dist, tau_prior){
        if(length(tau_prior) == 1L && is.na(tau_prior[1])) return(NA)
        list(dist = tau_dist, prior = tau_prior)
}

#-------------------------------------------------------------------------------
# Internal: null value of a contrast
#
# The single place that maps a contrast type to the value representing "no
# effect": ratio-type contrasts (mean_ratio, RR, OR / odds ratio,
# risk_ratio, rate_ratio, ...) test against 1, difference-type contrasts
# (mean_diff, pi_diff, risk_diff, rate_diff, ...) against 0. Keying on the
# substring "ratio" covers most future contrasts; the binomial "RR" (ratio
# of proportions) and "OR" (odds ratio) are named explicitly. Used by SCI plots,
# SIM_norm() and (future) the Poisson code.

.null_value <- function(contr){
        is_ratio <- length(contr) == 1L &&
                (grepl("ratio", contr, fixed = TRUE) || contr %in% c("RR", "OR"))
        if(is_ratio) 1 else 0
}

#-------------------------------------------------------------------------------
# Internal: human-readable label of a contrast type
#
# Singular form (for summary() headings) and plural / title-case form (for plot
# titles). Shared by summary.MultiMAP() and plot.MultiMAP() so normal and
# binomial endpoints are described consistently.

.contr_label <- function(contr){
        switch(as.character(contr),
               mean_diff  = "mean difference",
               mean_ratio = "mean ratio",
               pi_diff    = "difference of proportions",
               RR    = "ratio of proportions",
               OR         = "odds ratio",
               rate_diff  = "difference of rates",
               rate_ratio = "ratio of rates",
               "contrast")
}

.contr_title <- function(contr){
        switch(as.character(contr),
               mean_diff  = "Mean differences",
               mean_ratio = "Mean ratios",
               pi_diff    = "Differences of proportions",
               RR    = "Ratios of proportions",
               OR         = "Odds ratios",
               rate_diff  = "Differences of rates",
               rate_ratio = "Ratios of rates",
               "Contrasts")
}

#-------------------------------------------------------------------------------
# Internal helpers shared by the (multi-arm) borrowing pipeline
#-------------------------------------------------------------------------------

# Effective sample size of an RBesT mixture, robust to elir failures and capped
# at the total historical information 'cap'.
.ess_safe <- function(mix, cap){
        e <- try(round(ess(mix, method = "elir")), silent = TRUE)
        if(inherits(e, "try-error")) e <- round(ess(mix, method = "moment"))
        min(max(0, e), cap)
}

# Build a prediction interval (median/lower/upper) from an RBesT predictive
# mixture, honouring the one-sidedness 'limits' and dividing count quantiles by
# 'scale' (scale = 1 for the normal endpoint; n0 / offset0 for binom / pois so
# the interval is on the proportion / rate scale).
.pred_int_from <- function(pred_dist, limits, alpha, scale = 1){
        med <- unname(RBesT::qmix(pred_dist, 0.5)) / scale
        if(limits == "two.sided")
                c(median = med,
                  lower  = unname(RBesT::qmix(pred_dist, alpha / 2))     / scale,
                  upper  = unname(RBesT::qmix(pred_dist, 1 - alpha / 2)) / scale)
        else if(limits == "lower")
                c(median = med,
                  lower  = unname(RBesT::qmix(pred_dist, alpha)) / scale,
                  upper  = NA)
        else
                c(median = med,
                  lower  = NA,
                  upper  = unname(RBesT::qmix(pred_dist, 1 - alpha)) / scale)
}

# TRUE when histdat requests multi-arm borrowing (a named list of per-arm
# historical data sets) rather than the single (control-only) data frame.
.is_multi_arm <- function(histdat) is.list(histdat) && !is.data.frame(histdat)

# Validate a multi-arm 'histdat' (named list) against the current-trial groups.
.check_multi_histdat <- function(histdat, groups){
        if(is.null(names(histdat)) || any(names(histdat) == ""))
                stop("For multi-arm borrowing, 'histdat' must be a *named* list; ",
                     "each element's name is the arm it provides historical data for.")
        if(anyDuplicated(names(histdat)))
                stop("'histdat' has duplicated arm names.")
        miss <- setdiff(names(histdat), as.character(groups))
        if(length(miss))
                stop("These 'histdat' names do not match any arm in newdat's first column: ",
                     paste(miss, collapse = ", "), ".")
        if(!all(vapply(histdat, is.data.frame, logical(1))))
                stop("Every element of a multi-arm 'histdat' must be a data frame.")
        invisible(TRUE)
}

#-------------------------------------------------------------------------------
# Internal: frequentist, non-borrowing competitor (prior_est = "ignore_HCD").
#
# Analyses the current trial ALONE with a standard GLM + a single-step Dunnett-
# type multiple-comparison procedure (multcomp), using heteroscedastic
# (unequal-variance) per-arm variances -- matching how the Bayesian pipeline
# samples each arm with its own observed spread. Per-group means are taken on the
# working scale (identity / logit / log via a '~ 0 + group' fit), and the
# contrasts are formed with the same 'type'/'base' contrast matrix used by the
# borrowing analysis. Inference is large-sample (multivariate normal, df = 0),
# consistent with the Bayesian machinery treating the reference scale as known.
#
#   normal   -> mean differences (identity)
#   binomial -> odds ratios      (logit -> exp)
#   poisson  -> rate ratios      (log   -> exp)
#
# Returns the standardised 'comp' data frame (comp, estimate, lower.ci,
# upper.ci) and the contrast matrix. Fitting problems (e.g. glm errors on
# complete separation) propagate as errors and are handled by the caller.

.freq_competitor <- function(newdat, distr, type, base, alpha, alternative){

        # Group factor with the control (row 1) as the first level.
        grp      <- factor(as.character(newdat[[1]]), levels = as.character(newdat[[1]]))
        alt_glht <- switch(alternative, "lower" = "greater", "upper" = "less", "two.sided")

        if(distr == "normal"){
                theta     <- newdat$mean
                Sigma     <- diag(newdat$sd^2 / newdat$n, nrow = length(theta))
                nni       <- newdat$n
                exp_scale <- FALSE
        }
        else if(distr == "betabinomial"){
                dd    <- data.frame(ev = newdat[[2]], nonev = newdat[[3]], grp = grp)
                fit   <- stats::glm(cbind(ev, nonev) ~ 0 + grp, data = dd,
                                    family = stats::binomial())
                theta <- unname(fit$coefficients)
                Sigma <- stats::vcov(fit)
                nni   <- newdat[[2]] + newdat[[3]]
                exp_scale <- TRUE
        }
        else if(distr == "poisson"){
                dd    <- data.frame(ev = newdat[[2]], off = newdat[[3]], grp = grp)
                fit   <- stats::glm(ev ~ 0 + grp + offset(log(off)), data = dd,
                                    family = stats::poisson())
                theta <- unname(fit$coefficients)
                Sigma <- stats::vcov(fit)
                nni   <- newdat[[3]]
                exp_scale <- TRUE
        }

        names(theta)  <- levels(grp)
        names(nni)    <- levels(grp)

        K  <- multcomp::contrMat(nni, type = type, base = base)
        g  <- multcomp::glht(multcomp::parm(theta, Sigma), linfct = K,
                             alternative = alt_glht)
        ci <- stats::confint(g, level = 1 - alpha)$confint
        if(exp_scale) ci <- exp(ci)

        comp <- data.frame(row.names = seq_len(nrow(ci)),
                           comp     = rownames(ci),
                           estimate = unname(signif(ci[, 1], 3)),
                           lower.ci = unname(signif(ci[, 2], 3)),
                           upper.ci = unname(signif(ci[, 3], 3)))

        list(comp = comp, contr_mat = K)
}

#' Prior and prediction interval for a Poisson endpoint
#'
#' Builds robustified priors for the study arms of interest -- either meta-analytic-predictive (MAP) priors
#' (\code{prior_est = "MAP"}, via \code{RBesT::gMAP})
#' or empirical-Bayes gamma priors estimated from the historical data
#' (\code{prior_est = "EB"}) -- and derives the
#' prediction intervals for the rates in the study arms of interest.
#'
#' @param histdat external (historical) data for the study-arms of interest; 
#' borrowing for one study arm (e.g. control group): submit a data.frame with columns
#'   \code{events} (the counted endpoint) and \code{offset} (the exposure / number of
#'   experimental units). 
#'   To borrow for more than one arm, pass a \emph{named list}
#'   containing different data frames: each element supplies the external
#'    data for a particular study-arm. 
#' @param offset0 exposure (number of experimental units) in the study arms of interest
#'   for which the prediction intervals are derived. 
#'   Either a single value (used for every arm) or a named vector with one value per arm.
#' @param prior_est prior for the concurrent control rate: \code{"MAP"} (via
#'   \code{RBesT::gMAP}, \code{family = poisson}, fitted to a gamma mixture) or
#'   \code{"EB"} (an empirical-Bayes gamma prior following Tarone (1982)).
#' @param limits limits of the prediction interval: either \code{"two.sided"},
#'   \code{"lower"} or \code{"upper"}.
#' @param robust weight of the weakly-informative (robust) prior component added
#'   to robustify the prior.
#' @param tau_prior optional override for the \code{gMAP} between-study
#'   heterogeneity prior (\code{tau.prior} in \code{RBesT::gMAP}), on the
#'   log-rate scale. \code{NULL} (the default) uses a fixed
#'   \code{HalfNormal(0, 0.5)}. For a more diffuse / conservative
#'   prior (borrowing less under substantial heterogeneity) use \code{cbind(0, 1)},
#'   i.e. \code{HalfNormal(0, 1)}; supply any \code{cbind(0, scale)} to override.
#' @param tau_dist distribution family of the \code{gMAP} between-study
#'   heterogeneity prior (\code{tau.dist} in \code{RBesT::gMAP}); one of
#'   \code{"HalfNormal"} (the default), \code{"TruncNormal"}, \code{"Uniform"},
#'   \code{"Gamma"}, \code{"InvGamma"}, \code{"LogNormal"}, \code{"TruncCauchy"},
#'   \code{"Exp"} or \code{"Fixed"}. The built-in default \code{tau_prior}
#'   heuristic is calibrated for the half-normal, so when \code{tau_dist} is set
#'   to any other family a matching \code{tau_prior} must also be supplied.
#' @param beta_prior optional override for the \code{gMAP} intercept prior
#'   (\code{beta.prior} in \code{RBesT::gMAP}), on the log-rate scale. \code{NULL}
#'   (the default) uses \code{cbind(log(mean rate), 1)} -- a weakly-informative
#'   normal prior centred at the pooled historical log-rate (\code{mean rate} =
#'   mean of the per-study \code{events / offset}). For a more diffuse prior use an
#'   SD of 2, e.g. \code{cbind(log(mean rate), 2)}.
#' @param alpha 1 - coverage level of the prediction interval.
#' @param n_mcmc total number of Monte-Carlo samples to draw.
#' @param seed optional integer; if supplied, \code{set.seed(seed)} is called so
#'   the \code{gMAP} MCMC (for \code{prior_est = "MAP"}) is reproducible.
#'   \code{NULL} (default) leaves the RNG state untouched.
#' @param all_zero length-2 numeric giving the fallback gamma shape and rate
#'   \code{c(a, b)} used for the empirical-Bayes prior when all historical counts
#'   are zero.
#'
#' @return A MultiMAP object: a named list with S3 class
#'   \code{c("PI", "MultiMAP_pois", "MultiMAP")} and the following slots:
#'   \describe{
#'     \item{\code{comp}}{\code{NA} (only \code{\link{SCI_pois}} fills this).}
#'     \item{\code{contr}}{\code{NA} (only \code{\link{SCI_pois}} fills this).}
#'     \item{\code{histdat}}{the historical control data used.}
#'     \item{\code{newdat}}{\code{NA} (no current trial groups in a
#'       prediction-interval fit).}
#'     \item{\code{distr}}{the endpoint family, \code{"poisson"}.}
#'     \item{\code{prior_par}}{\code{NA} (only \code{\link{SCI_pois}} fills this).}
#'     \item{\code{weights}}{data frame of the prior weights of the robustified
#'       mixture's components (posterior weights are \code{NA}).}
#'     \item{\code{eff_n_map}}{effective sample size of the (non-robustified)
#'       prior.}
#'     \item{\code{eff_n}}{effective sample size of the robustified prior.}
#'     \item{\code{eff_n_post}}{\code{NA} (no posterior is formed).}
#'     \item{\code{estimates}}{\code{NA} (only \code{\link{SCI_pois}} fills this).}
#'     \item{\code{contr_mat}}{\code{NA} (only \code{\link{SCI_pois}} fills this).}
#'     \item{\code{prior}}{the (non-robustified) prior (RBesT gamma mixture).}
#'     \item{\code{prior_rob}}{the robustified prior (RBesT gamma mixture).}
#'     \item{\code{posterior}}{\code{NA} (no posterior is formed).}
#'     \item{\code{pred_dist}}{the prior-predictive distribution for the number of
#'       events in a future concurrent control group, an RBesT mixture object
#'       built with \code{\link[RBesT]{preddist}}.}
#'     \item{\code{pred_int}}{the prediction interval for the concurrent control
#'       rate (named \code{median}, \code{lower}, \code{upper}).}
#'     \item{\code{shrinkage}}{\code{NA} (only \code{\link{SCI_pois}} fills this).}
#'     \item{\code{simpost}}{\code{NA} (only \code{\link{SCI_pois}} fills this).}
#'     \item{\code{tau_prior}}{a list with the gMAP between-study heterogeneity
#'       prior distribution name (\code{dist}) and its parameters (\code{prior})
#'       actually used (\code{NA} for the EB prior).}
#'     \item{\code{beta_prior}}{the gMAP intercept prior actually used
#'       (\code{NA} for the EB prior).}
#'     \item{\code{alpha}}{the \code{alpha} used.}
#'     \item{\code{seed}}{the \code{seed} used (\code{NA} if none).}
#'   }
#'   Use \code{\link{summary.MultiMAP}}, \code{\link{get_output}} and
#'   \code{\link{report_pois}} on it.
#'
#' @seealso \code{\link{SCI_pois}}, \code{\link{report_pois}}
#'
#' @references
#' EFSA (European Food Safety Authority) (2025). Use and reporting of historical
#' control data for regulatory studies. \emph{EFSA Journal}.
#' \doi{10.2903/j.efsa.2025.9576}
#'
#' Menssen, M., Kneuer, C., Akyianu, G., Roever, C., Friede, T., &
#' Schaarschmidt, F. (2026). Including historical control data in simultaneous
#' inference for pre-clinical multi-arm studies. \emph{arXiv preprint}
#' arXiv:2603.11730. \doi{10.48550/arXiv.2603.11730}
#'
#' Schmidli, H., Gsteiger, S., Roychoudhury, S., O'Hagan, A., Spiegelhalter, D.,
#' & Neuenschwander, B. (2014). Robust meta-analytic-predictive priors in
#' clinical trials with historical control information. \emph{Biometrics,
#' 70}(4), 1023--1032. \doi{10.1111/biom.12242}
#'
#' Roever, C., Bender, R., Dias, S., Schmid, C. H., Schmidli, H., Sturtz, S.,
#' Weber, S., & Friede, T. (2021). On weakly informative prior distributions for
#' the heterogeneity parameter in Bayesian random-effects meta-analysis.
#' \emph{Research Synthesis Methods, 12}(4), 448--474. \doi{10.1002/jrsm.1475}
#'
#' Tarone, R. E. (1982). The use of historical control information in testing for
#' a trend in Poisson means. \emph{Biometrics, 38}(2), 457--462.
#' \doi{10.2307/2530459}
#'
#' @examples
#' \donttest{
#' # Egg-laying count data from EFSA (2025), Annex C: number of eggs laid by a
#' # control group of 12 hens over 10 weeks. HCD = 106 historical control studies
#' # summarised by their mean eggs/hen (Table 3). Poisson mapping: offset = 12
#' # hens, events = total eggs = round(mean eggs/hen * 12), so rate = eggs/hen.
#' hcd_mean <- c(
#'   43.1, 40.5, 48.2, 30.0, 44.9, 44.8, 42.0, 39.2, 45.9, 48.4, 41.0, 48.3,
#'   39.4, 34.9, 36.6, 52.2, 36.5, 39.2, 37.2, 38.6, 43.6, 44.8, 46.9, 51.2,
#'   38.2, 33.1, 43.6, 46.2, 49.8, 35.9, 37.2, 35.9, 40.8, 55.5, 47.8, 47.8,
#'   44.1, 44.3, 52.0, 28.2, 30.2, 41.4, 40.1, 32.3, 43.5, 41.3, 47.5, 39.5,
#'   43.4, 44.2, 44.3, 36.7, 48.4, 39.2, 40.8, 35.1, 30.8, 47.0, 38.4, 43.5,
#'   39.2, 39.4, 43.5, 45.6, 40.2, 43.4, 37.8, 24.7, 29.4, 32.2, 43.4, 41.6,
#'   35.2, 44.2, 38.7, 38.8, 42.8, 37.8, 35.7, 38.6, 38.7, 41.3, 26.8, 41.0,
#'   43.5, 37.2, 39.6, 43.5, 33.8, 31.7, 41.4, 29.5, 40.8, 23.9, 40.6, 42.8,
#'   33.9, 35.8, 34.4, 35.7, 31.7, 45.1, 30.8, 36.5, 51.3, 54.0)
#' histdat <- data.frame(events = round(hcd_mean * 12), offset = 12)
#'
#' # Prediction interval for a future concurrent control group of 12 hens
#' m_pi <- PI_pois(histdat = histdat, offset0 = 12, prior_est = "MAP",
#'                 seed = 81771)
#' summary(m_pi)
#' }
#' @export
PI_pois <- function(histdat,
                    offset0,
                    prior_est = "MAP",
                    limits = "two.sided",
                    robust = 0.2,
                    tau_prior = NULL,
                    tau_dist = "HalfNormal",
                    beta_prior = NULL,
                    alpha = 0.05,
                    n_mcmc = 50000,
                    seed = NULL,
                    all_zero = c(0.0001, 1)){

        # Optional reproducibility: seeds gMAP()'s MCMC (MAP path).
        if(!is.null(seed)) set.seed(seed)

        # User-facing one-sidedness vocabulary: two.sided / lower / upper.
        limits <- match.arg(limits, c("two.sided", "lower", "upper"))

        # gMAP heterogeneity-prior distribution (tau.dist); non-half-normal
        # families need an explicit tau_prior.
        tau_dist <- match.arg(tau_dist, .tau_dist_choices)
        .check_tau_dist(tau_dist, tau_prior)

        # Multi-arm: 'histdat' a named list gives one prediction interval per arm.
        multi_arm <- .is_multi_arm(histdat)

        #-----------------------------------------------------------------------
        # Input checks

        if(!multi_arm && (!is.data.frame(histdat) || ncol(histdat) < 2))
                stop("'histdat' must be a data frame with at least 2 columns (events, offset).")
        if(!is.numeric(offset0) || !all(offset0 > 0))
                stop("'offset0' must be positive (a single value or one per borrowed arm).")
        if(!multi_arm && length(offset0) != 1)
                stop("'offset0' must be a single positive number.")
        if(!prior_est %in% c("MAP", "EB"))
                stop("'prior_est' must be either 'MAP' or 'EB'.")
        if(!is.numeric(robust) || length(robust) != 1 || robust < 0 || robust > 1)
                stop("'robust' must be a single number in [0, 1].")
        if(!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1)
                stop("'alpha' must be a single number in (0, 1).")

        #-----------------------------------------------------------------------
        # Multi-arm: one MAP/EB prior + prediction interval per named arm.
        if(multi_arm){
                if(is.null(names(histdat)) || any(names(histdat) == "") ||
                   anyDuplicated(names(histdat)) ||
                   !all(vapply(histdat, is.data.frame, logical(1))))
                        stop("For multi-arm borrowing, 'histdat' must be a named list of data frames.")
                return(.pi_pois_multiarm(histdat, offset0, prior_est, robust,
                                         tau_prior, tau_dist, beta_prior, limits, alpha, n_mcmc,
                                         all_zero, seed))
        }

        # Effective gMAP hyperpriors actually used (filled in the MAP block; NA
        # for the EB path). Stored in the output.
        tau_p  <- NA
        beta_p <- NA

        #=======================================================================
        # Prior specification -- same as SCI_pois(): MAP (gMAP) or EB (Tarone)
        #=======================================================================

        if(prior_est == "MAP"){

                histdat_map <- cbind(data.frame(study = factor(seq_len(nrow(histdat)))),
                                     histdat)
                colnames(histdat_map)[2:3] <- c("events", "offset")

                # Pooled historical (grand) rate; its log centres the intercept prior.
                lambda_bar <- mean(histdat_map$events / histdat_map$offset)

                if(!is.finite(lambda_bar) || lambda_bar <= 0)
                        stop("The MAP prior needs positive historical counts; use prior_est = 'EB' for all-zero data.")

                # Default gMAP hyperpriors on the log-rate scale, kept identical to
                # SCI_pois() so PI_pois() and SCI_pois() share the same MAP prior:
                #  - tau ~ HalfNormal(0, 0.5): a fixed, weakly-informative between-study
                #    heterogeneity prior (heterogeneity is a coefficient of variation on
                #    the log-rate scale, ~ invariant to the rate; Roever et al. 2021).
                #    cbind(0, 1) gives a more diffuse alternative.
                #  - intercept ~ Normal(log(lambda_bar), 1): centred at the pooled
                #    historical log-rate with a weakly-informative SD (use 2 for a more
                #    diffuse prior). Both can be overridden via tau_prior / beta_prior.
                tau_p  <- if(is.null(tau_prior))  cbind(0, 0.5)             else tau_prior
                beta_p <- if(is.null(beta_prior)) cbind(log(lambda_bar), 1) else beta_prior

                map_mcmc <- gMAP(events ~ 1 + offset(log(offset)) | study,
                                 data       = histdat_map,
                                 family     = poisson,
                                 tau.dist   = tau_dist,
                                 tau.prior  = tau_p,
                                 beta.prior = beta_p)

                map <- automixfit(map_mcmc, type = "gamma")
        }

        if(prior_est == "EB"){

                if(all(histdat[, 1] == 0)){
                        a_est <- all_zero[1]
                        b_est <- all_zero[2]
                }
                else{
                        # Tarone (1982) gamma prior Gamma(p, p / mu); shape p by the
                        # method of moments, capped at the pooled data for
                        # equi-/under-dispersed or very homogeneous counts. See
                        # SCI_pois() for the full rationale.
                        y_h    <- histdat[, 1]
                        off_h  <- histdat[, 2]
                        lambda_hat <- sum(y_h) / sum(off_h)
                        mu_h       <- lambda_hat * off_h
                        denom      <- sum((y_h - mu_h)^2 - mu_h)

                        if(denom <= 0){
                                a_est <- sum(y_h)
                                b_est <- sum(off_h)
                                warning("HCD is equi-/under-dispersed: the prior was capped at the pooled historical data (full pooling).")
                        }
                        else{
                                p_hat <- sum(mu_h^2) / denom
                                a_est <- p_hat
                                b_est <- p_hat / lambda_hat
                                if(b_est > sum(off_h)){
                                        a_est <- sum(y_h)
                                        b_est <- sum(off_h)
                                }
                        }
                }

                map <- mixgamma(c(1, a_est, b_est), param = "ab", likelihood = "poisson")
        }

        #=======================================================================
        # Robustify + prediction interval
        #=======================================================================

        prior_mean <- unname(summary(map)["mean"])
        map_robust <- robustify(map, weight = robust, mean = prior_mean)

        #-----------------------------------------------------------------------
        # Prediction interval for the concurrent control rate.
        # preddist() on the (non-robustified) prior gives the gamma-Poisson
        # predictive for the number of events over the future exposure; dividing
        # the count quantiles by that exposure puts the interval on the rate scale.

        pred_dist <- RBesT::preddist(map, n = offset0)

        if(limits == "two.sided"){
                pred_int <- c(median = unname(RBesT::qmix(pred_dist, 0.5))        / offset0,
                              lower  = unname(RBesT::qmix(pred_dist, alpha/2))     / offset0,
                              upper  = unname(RBesT::qmix(pred_dist, 1 - alpha/2)) / offset0)
        }

        if(limits == "lower"){
                pred_int <- c(median = unname(RBesT::qmix(pred_dist, 0.5))   / offset0,
                              lower  = unname(RBesT::qmix(pred_dist, alpha))  / offset0,
                              upper  = NA)
        }

        if(limits == "upper"){
                pred_int <- c(median = unname(RBesT::qmix(pred_dist, 0.5))        / offset0,
                              lower  = NA,
                              upper  = unname(RBesT::qmix(pred_dist, 1 - alpha))  / offset0)
        }

        #-----------------------------------------------------------------------
        # Effective sample size of the (non-robustified) prior

        eff_n_map <- try(round(ess(map, method = "elir")), silent = TRUE)

        if(inherits(eff_n_map, "try-error")){

                eff_n_map <- round(ess(map, method = "moment"))

                warning("Estimation of ESS via elir-method does not work.\n ESS was estimated based on method of moments")
        }

        eff_n_map <- min(max(0, eff_n_map), sum(histdat[, 2]))

        #-----------------------------------------------------------------------
        # Effective sample size of the robustified prior

        eff_sample_size <- try(round(ess(map_robust, method = "elir")), silent = TRUE)

        if(inherits(eff_sample_size, "try-error")){

                eff_sample_size <- round(ess(map_robust, method = "moment"))

                warning("Estimation of ESS via elir-method does not work.\n ESS was estimated based on method of moments")
        }

        eff_sample_size <- min(max(0, eff_sample_size), sum(histdat[, 2]))

        #-----------------------------------------------------------------------
        weights <- data.frame(w_prior = c(1 - robust, robust),
                              w_post  = c(NA, NA))

        #-----------------------------------------------------------------------
        # Output

        out_list <- list("comp"       = NA,
                         "contr"      = NA,
                         "histdat"    = histdat,
                         "newdat"     = NA,
                         "distr"      = "poisson",
                         "prior_par"  = NA,
                         "weights"    = weights,
                         "eff_n_map"  = eff_n_map,
                         "eff_n"      = eff_sample_size,
                         "eff_n_post" = NA,
                         "estimates"  = NA,
                         "contr_mat"  = NA,
                         "prior"      = map,
                         "prior_rob"  = map_robust,
                         "posterior"  = NA,
                         "pred_dist"  = pred_dist,
                         "pred_int"   = pred_int,
                         "shrinkage"  = NA,
                         "simpost"    = NA,
                         "tau_prior"  = .tau_prior_out(tau_dist, tau_p),
                         "beta_prior" = beta_p,
                         "alpha"      = alpha,
                         "seed"       = if(is.null(seed)) NA_integer_ else seed)

        out_s3 <- structure(out_list,
                            class = c("PI", "MultiMAP_pois", "MultiMAP"))

        return(out_s3)
}


#-------------------------------------------------------------------------------
# Internal: MAP/EB prior + prediction interval for ONE arm (Poisson, no
# posterior update). Reuses .map_prior_pois() (defined in SCI_pois.R).
#-------------------------------------------------------------------------------
.pi_pois_arm <- function(hd, offset0, prior_est, robust, tau_prior, tau_dist, beta_prior,
                         limits, alpha, all_zero){

        pr  <- .map_prior_pois(hd, prior_est, tau_prior, tau_dist, beta_prior, all_zero)
        map <- pr$map
        prior_mean <- unname(summary(map)["mean"])
        map_robust <- robustify(map, weight = robust, mean = prior_mean)
        pred_dist  <- RBesT::preddist(map, n = offset0)
        pred_int   <- .pred_int_from(pred_dist, limits, alpha, scale = offset0)
        cap        <- sum(hd[, 2])

        list(prior      = map,
             prior_rob  = map_robust,
             pred_dist  = pred_dist,
             pred_int   = pred_int,
             eff_n_map  = .ess_safe(map,        cap),
             eff_n      = .ess_safe(map_robust, cap),
             weights    = data.frame(w_prior = c(1 - robust, robust), w_post = c(NA, NA)),
             tau_prior  = .tau_prior_out(tau_dist, pr$tau_p),
             beta_prior = pr$beta_p)
}

#-------------------------------------------------------------------------------
# Internal: multi-arm PI_pois(). One prediction interval per named arm.
#-------------------------------------------------------------------------------
.pi_pois_multiarm <- function(histdat, offset0, prior_est, robust,
                              tau_prior, tau_dist, beta_prior, limits, alpha, n_mcmc,
                              all_zero, seed){

        if(!is.null(seed)) set.seed(seed)
        arms <- names(histdat)
        offv <- if(length(offset0) == 1) stats::setNames(rep(offset0, length(arms)), arms) else offset0[arms]
        if(any(is.na(offv)))
                stop("A named 'offset0' must supply a value for every arm in 'histdat'.")

        b <- lapply(arms, function(a)
                    .pi_pois_arm(histdat[[a]], offv[[a]], prior_est, robust,
                                 tau_prior, tau_dist, beta_prior, limits, alpha, all_zero))
        names(b) <- arms
        pick  <- function(nm) lapply(b, `[[`, nm)
        pickv <- function(nm) vapply(b, `[[`, numeric(1), nm)
        na_v  <- stats::setNames(rep(NA_real_, length(arms)), arms)

        out_list <- list("comp"        = NA, "contr" = NA, "histdat" = histdat,
                         "newdat"      = NA, "distr" = "poisson",
                         "borrow_arms" = arms, "prior_par" = NA,
                         "weights"     = pick("weights"),
                         "eff_n_map"   = pickv("eff_n_map"),
                         "eff_n"       = pickv("eff_n"),
                         "eff_n_post"  = na_v,
                         "estimates"   = NA, "contr_mat" = NA,
                         "prior"       = pick("prior"),
                         "prior_rob"   = pick("prior_rob"),
                         "posterior"   = NA,
                         "pred_dist"   = pick("pred_dist"),
                         "pred_int"    = pick("pred_int"),
                         "shrinkage"   = NA, "simpost" = NA,
                         "tau_prior"   = pick("tau_prior"),
                         "beta_prior"  = pick("beta_prior"),
                         "alpha"       = alpha,
                         "seed"        = if(is.null(seed)) NA_integer_ else seed)

        structure(out_list, class = c("PI", "MultiMAP_pois", "MultiMAP"))
}

#' MAP prior and prediction interval for a normal endpoint
#'
#' Builds robustified meta-analytic-predictive (MAP) priors (via \code{RBesT::gMAP}) 
#' and derives the prediction intervals for the means of the study arms of interest.
#'
#' @param histdat external (historical) data for the study-arms of interest; 
#' borrowing for one study arm (e.g. control group): data.frame with one row per study, 
#' columns \code{mean}, \code{sd}, \code{n}. To borrow for more than one arm, 
#' pass a \emph{named list} containing different data frames: each element supplies 
#' the external data for a particular study-arm. 
#' @param n0 planned sample size of the (future) concurrent control group
#'   (\eqn{m = 0}) for which the prediction interval is derived. For a named-list
#'   \code{histdat}, either a single value (used for every arm) or a named vector
#'   with one value per arm.
#' @param prior_est currently, only \code{"MAP"} is supported: a
#'   meta-analytic-predictive prior via \code{RBesT::gMAP} (\code{family =
#'   gaussian}), fitted to a normal mixture.
#' @param limits limits of the prediction interval: either \code{"two.sided"},
#'   \code{"lower"} or \code{"upper"}.
#' @param robust weight of the weakly-informative (robust) prior component added
#'   to robustify the MAP prior.
#' @param sigma reference scale for the RBesT Gaussian machinery. \code{NULL}
#'   (the default) estimates it from the data (n-weighted pooled sd of the
#'   historical controls); a supplied value is used as a fixed known reference
#'   scale throughout.
#' @param tau_prior optional override for the \code{gMAP} between-study
#'   heterogeneity prior (\code{tau.prior} in \code{RBesT::gMAP}). \code{NULL}
#'   (the default) uses the built-in \code{HalfNormal(0, sigma/2)}; supply e.g.
#'   \code{cbind(0, scale)}. A larger scale widens the MAP prior and its
#'   prediction interval -- useful when few historical studies leave tau weakly
#'   identified.
#' @param tau_dist distribution family of the \code{gMAP} between-study
#'   heterogeneity prior (\code{tau.dist} in \code{RBesT::gMAP}); one of
#'   \code{"HalfNormal"} (the default), \code{"TruncNormal"}, \code{"Uniform"},
#'   \code{"Gamma"}, \code{"InvGamma"}, \code{"LogNormal"}, \code{"TruncCauchy"},
#'   \code{"Exp"} or \code{"Fixed"}. The built-in default \code{tau_prior}
#'   heuristic is calibrated for the half-normal, so when \code{tau_dist} is set
#'   to any other family a matching \code{tau_prior} must also be supplied.
#' @param beta_prior optional override for the \code{gMAP} intercept prior
#'   (\code{beta.prior} in \code{RBesT::gMAP}). \code{NULL} (the default) uses
#'   \code{cbind(grand mean, sigma)}.
#' @param alpha 1 - coverage level of the prediction interval.
#' @param n_mcmc total number of Monte-Carlo samples to draw.
#' @param seed optional integer; if supplied, \code{set.seed(seed)} is called so
#'   the \code{gMAP} MCMC is reproducible. \code{NULL} (default) leaves the RNG
#'   state untouched.
#'
#' @return A MultiMAP object: a named list with S3 class
#'   \code{c("PI", "MultiMAP_norm", "MultiMAP")} and the following slots:
#'   \describe{
#'     \item{\code{comp}}{\code{NA} (only \code{\link{SCI_norm}} fills this).}
#'     \item{\code{contr}}{\code{NA} (only \code{\link{SCI_norm}} fills this).}
#'     \item{\code{histdat}}{the historical control data used (mean, sd, n).}
#'     \item{\code{newdat}}{\code{NA} (no current trial groups in a
#'       prediction-interval fit).}
#'     \item{\code{distr}}{the endpoint family, \code{"normal"}.}
#'     \item{\code{prior_par}}{\code{NA} (only \code{\link{SCI_norm}} fills this).}
#'     \item{\code{weights}}{data frame of the prior weights of the robustified
#'       mixture's components (posterior weights are \code{NA}).}
#'     \item{\code{eff_n_map}}{effective sample size of the (non-robustified) MAP
#'       prior.}
#'     \item{\code{eff_n}}{effective sample size of the robustified MAP prior.}
#'     \item{\code{eff_n_post}}{\code{NA} (no posterior is formed).}
#'     \item{\code{estimates}}{\code{NA} (only \code{\link{SCI_norm}} fills this).}
#'     \item{\code{contr_mat}}{\code{NA} (only \code{\link{SCI_norm}} fills this).}
#'     \item{\code{prior}}{the MAP prior (RBesT normal mixture).}
#'     \item{\code{prior_rob}}{the robustified MAP prior (RBesT normal mixture).}
#'     \item{\code{posterior}}{\code{NA} (no posterior is formed).}
#'     \item{\code{pred_dist}}{the MAP-based predictive distribution for a future
#'       concurrent control mean, an RBesT mixture object built with
#'       \code{\link[RBesT]{preddist}}.}
#'     \item{\code{pred_int}}{the prediction interval (named \code{median},
#'       \code{lower}, \code{upper}).}
#'     \item{\code{shrinkage}}{\code{NA} (only \code{\link{SCI_norm}} fills this).}
#'     \item{\code{simpost}}{\code{NA} (only \code{\link{SCI_norm}} fills this).}
#'     \item{\code{tau_prior}}{a list with the gMAP between-study heterogeneity
#'       prior distribution name (\code{dist}) and its parameters (\code{prior})
#'       actually used.}
#'     \item{\code{beta_prior}}{the gMAP intercept prior actually used.}
#'     \item{\code{alpha}}{the \code{alpha} used.}
#'     \item{\code{seed}}{the \code{seed} used (\code{NA} if none).}
#'   }
#'   Use \code{\link{plot.MultiMAP}}, \code{\link{summary.MultiMAP}},
#'   \code{\link{get_output}} and \code{\link{report_norm}} on it.
#'
#' @seealso \code{\link{SCI_norm}}, \code{\link{report_norm}}
#'
#' @references
#' EFSA (European Food Safety Authority) (2025). Use and reporting of historical
#' control data for regulatory studies. \emph{EFSA Journal}.
#' \doi{10.2903/j.efsa.2025.9576}
#'
#' Schmidli, H., Gsteiger, S., Roychoudhury, S., O'Hagan, A., Spiegelhalter, D.,
#' & Neuenschwander, B. (2014). Robust meta-analytic-predictive priors in
#' clinical trials with historical control information. \emph{Biometrics,
#' 70}(4), 1023--1032. \doi{10.1111/biom.12242}
#'
#' @examples
#' \donttest{
#' # Historical control data (diet route) for a clinical biochemical endpoint
#' # (mmol/l) in male rats, from EFSA (2025): study sizes from Table 13, per-study
#' # means and standard errors from Figure 4 (the 26 diet-route studies).
#' hcd <- data.frame(
#'   mean = c(1.11, 1.00, 1.37, 1.19, 1.34, 1.48, 1.51, 1.08, 1.44, 1.23, 1.38,
#'            1.44, 1.10, 1.77, 1.52, 1.22, 1.34, 1.37, 1.31, 1.24, 1.42, 1.35,
#'            1.14, 1.28, 1.35, 1.46),
#'   se   = c(0.143, 0.134, 0.164, 0.149, 0.161, 0.173, 0.176, 0.140, 0.169,
#'            0.152, 0.165, 0.169, 0.142, 0.198, 0.177, 0.151, 0.161, 0.164,
#'            0.159, 0.153, 0.168, 0.162, 0.145, 0.156, 0.162, 0.171),
#'   n    = c(5, 5, 5, 5, 10, 10, 5, 5, 10, 10, 10, 5, 5, 5, 5, 5, 5, 5, 5, 5,
#'            5, 5, 5, 5, 5, 5))
#'
#' # SCI_norm()/PI_norm() need the per-study SD: sd = se * sqrt(n)
#' histdat <- data.frame(mean = hcd$mean, sd = hcd$se * sqrt(hcd$n), n = hcd$n)
#'
#' # Prediction interval for a future concurrent control group of 5 animals
#' m_pi <- PI_norm(histdat = histdat, n0 = 5, seed = 81771)
#' summary(m_pi)
#' }
#' @export
PI_norm <- function(histdat,
                    n0,
                    prior_est = "MAP",
                    limits = "two.sided",
                    robust = 0.2,
                    sigma = NULL,
                    tau_prior = NULL,
                    tau_dist = "HalfNormal",
                    beta_prior = NULL,
                    alpha = 0.05,
                    n_mcmc = 50000,
                    seed = NULL){

        # Optional reproducibility: seeds gMAP()'s MCMC.
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

        if(!multi_arm && (!is.data.frame(histdat) || ncol(histdat) < 3))
                stop("'histdat' must be a data frame with at least 3 columns (mean, sd, n).")
        if(!is.numeric(n0) || !all(n0 > 0))
                stop("'n0' must be positive (a single value or one per borrowed arm).")
        if(!multi_arm && length(n0) != 1)
                stop("'n0' must be a single positive number.")
        if(!is.numeric(robust) || length(robust) != 1 || robust < 0 || robust > 1)
                stop("'robust' must be a single number in [0, 1].")
        if(!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1)
                stop("'alpha' must be a single number in (0, 1).")

        # Only the MAP prior is supported (the pooled approach has been removed).
        if(prior_est != "MAP"){
                stop("PI_norm() supports only prior_est = 'MAP'.")
        }

        #-----------------------------------------------------------------------
        # Multi-arm: one MAP prior + prediction interval per named arm.
        if(multi_arm){
                if(is.null(names(histdat)) || any(names(histdat) == "") ||
                   anyDuplicated(names(histdat)) ||
                   !all(vapply(histdat, is.data.frame, logical(1))))
                        stop("For multi-arm borrowing, 'histdat' must be a named list of data frames.")
                return(.pi_norm_multiarm(histdat, n0, robust, sigma,
                                         tau_prior, tau_dist, beta_prior, limits, alpha, seed))
        }

        #-----------------------------------------------------------------------
        # Column bookkeeping (mean, sd, n)

        colnames(histdat)[1:3] <- c("mean", "sd", "n")
        
        #-----------------------------------------------------------------------
        # Reference scale (sigma) for the RBesT Gaussian machinery
        
        # Unit-information standard deviation (UISD) of the historical controls:
        # the single observation-scale reference SD implied by the per-study
        # sample sizes and SDs. Same estimator as bayesmeta::uisd() and the
        # reference scale RBesT's gaussian gMAP derives internally, i.e. an
        # n-weighted harmonic pooling of the study variances
        #   sigma_u = sqrt( sum(n_h) / sum(n_h / sd_h^2) ).
        sigma_pool <- sqrt(sum(histdat$n) /
                                   sum(histdat$n / histdat$sd^2))
        
        if(!is.finite(sigma_pool) | sigma_pool <= 0){
                sigma_pool <- mean(histdat$sd)
        }
        
        # Reference scale: user-supplied value if given, else the pooled sd.
        # In the MAP path with sigma = NULL this is later refined to the
        # reference scale that gMAP estimates from the data.
        sigma_ref <- if(is.null(sigma)) sigma_pool else sigma
        
        # Robust prior mean: grand (n-weighted) mean of the historical controls
        robust_mean <- sum(histdat$n * histdat$mean) / sum(histdat$n)

        # Effective gMAP hyperpriors actually used (filled in the MAP block).
        tau_p  <- NA
        beta_p <- NA

        #-----------------------------------------------------------------------

        if(prior_est == "MAP"){

                histdat_map <- cbind(data.frame(study = factor(seq_len(nrow(histdat)))),
                                     histdat)


                # Per-study standard error of the mean
                histdat_map$se <- histdat_map$sd / sqrt(histdat_map$n)
                
                #---------------------------------------------------------------
                # MAP prior via gMAP (gaussian family).
                # For a normal endpoint gMAP expects the response as a two-column
                # matrix cbind(mean, se); 'weights = n' supplies the per-study unit
                # counts used to estimate the reference scale sigma.
                # tau ~ HalfNormal(0, sigma/2) follows the n_infinity heuristic.
                # beta.prior is the RBesT-recommended unit-information prior on the
                # intercept: SD = sigma (one observation's worth of information).
                # Unlike the RBesT normal vignette -- whose response is a
                # change-from-baseline effect and is therefore centred at 0 -- our
                # response is the raw group mean, so the prior is centred at the
                # historical grand mean rather than at 0.
                
                tau_p  <- if(is.null(tau_prior))  cbind(0, sigma_ref / 2)        else tau_prior
                beta_p <- if(is.null(beta_prior)) cbind(robust_mean, sigma_ref)  else beta_prior

                map_mcmc <- gMAP(cbind(mean, se) ~ 1 | study,
                                 data       = histdat_map,
                                 family     = gaussian,
                                 weights    = n,
                                 tau.dist   = tau_dist,
                                 tau.prior  = tau_p,
                                 beta.prior = beta_p)
                
                #---------------------------------------------------------------
                # Derive normal mixture
                
                map <- automixfit(map_mcmc, type = "norm")
                
                # Reference scale: enforce a user-supplied sigma, otherwise adopt
                # the data-driven scale that gMAP/automixfit estimated.
                if(!is.null(sigma)){
                        sigma(map) <- sigma_ref
                }
                else{
                        sigma_ref <- sigma(map)
                }
        }
        
        #-----------------------------------------------------------------------
        # Derive prediction interval from MAP
        
        pred_dist <- RBesT::preddist(map, n = n0, sigma=sigma_ref)     # uses sigma(mapmix) as reference scale
        
        if(limits == "two.sided"){
                pred_int <- c(median = unname(RBesT::qmix(pred_dist, 0.5)),
                              lower  = unname(RBesT::qmix(pred_dist, alpha/2)),
                              upper  = unname(RBesT::qmix(pred_dist, 1 - alpha/2)))
                
        }
        
        if(limits == "lower"){
                pred_int <- c(median = unname(RBesT::qmix(pred_dist, 0.5)),
                              lower  = unname(RBesT::qmix(pred_dist, alpha)),
                              upper  = NA)
                
        }
        
        if(limits == "upper"){
                pred_int <- c(median = unname(RBesT::qmix(pred_dist, 0.5)),
                              lower  = NA,
                              upper  = unname(RBesT::qmix(pred_dist, 1 - alpha)))

        }
        
        
        #-----------------------------------------------------------------------
        # Add weakly-informative (robust) component.
        # Following the RBesT normal vignette, the robust component is centred on
        # the mean of the (informative) MAP prior itself -- not on the raw grand
        # mean -- so that robustification only adds a heavy, vague tail without
        # shifting the prior's location.
        
        map_mean <- summary(map)["mean"]
        
        map_robust <- robustify(map,
                                weight = robust,
                                mean   = map_mean,
                                sigma  = sigma_ref)
        
        map_robust_df <- data.frame(cbind(t(map_robust)))
        
        #-----------------------------------------------------------------------
        # Effective sample size of the (non-robustified) MAP prior

        eff_n_map <- try(round(ess(map, method = "elir")), silent = TRUE)

        if(inherits(eff_n_map, "try-error")){

                eff_n_map <- round(ess(map, method = "moment"))

                warning("Estimation of ESS via elir-method does not work.\n ESS was estimated based on method of moments")
        }

        eff_n_map <- min(max(0, eff_n_map),
                         sum(histdat$n))

        #-----------------------------------------------------------------------
        # Effective sample size of the robustified MAP prior

        eff_sample_size <- try(round(ess(map_robust, method = "elir")), silent = TRUE)

        if(inherits(eff_sample_size, "try-error")){

                eff_sample_size <- round(ess(map_robust, method = "moment"))

                warning("Estimation of ESS via elir-method does not work.\n ESS was estimated based on method of moments")
        }

        eff_sample_size <- min(max(0, eff_sample_size),
                               sum(histdat$n))
        
        
        #-----------------------------------------------------------------------
        weights <- data.frame(w_prior = c(1 - robust, robust),
                              w_post  = c(NA, NA))
        
        #-----------------------------------------------------------------------
        # Output
        
        out_list <- list("comp"       = NA,
                         "contr"      = NA,
                         "histdat"    = histdat,
                         "newdat"     = NA,
                         "distr"      = "normal",
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
                            class = c("PI", "MultiMAP_norm", "MultiMAP"))

        return(out_s3)
}


#-------------------------------------------------------------------------------
# Internal: MAP prior + prediction interval for ONE arm (normal endpoint, no
# posterior update). Mirrors the single-arm PI_norm() prior block.
#-------------------------------------------------------------------------------
.pi_norm_arm <- function(hd, n0, robust, sigma, tau_prior, tau_dist, beta_prior, limits, alpha){

        colnames(hd)[1:3] <- c("mean", "sd", "n")

        sigma_pool <- sqrt(sum(hd$n) / sum(hd$n / hd$sd^2))
        if(!is.finite(sigma_pool) || sigma_pool <= 0) sigma_pool <- mean(hd$sd)
        sigma_ref   <- if(is.null(sigma)) sigma_pool else sigma
        robust_mean <- sum(hd$n * hd$mean) / sum(hd$n)

        hd_map    <- cbind(data.frame(study = factor(seq_len(nrow(hd)))), hd)
        hd_map$se <- hd_map$sd / sqrt(hd_map$n)

        tau_p  <- if(is.null(tau_prior))  cbind(0, sigma_ref / 2)       else tau_prior
        beta_p <- if(is.null(beta_prior)) cbind(robust_mean, sigma_ref) else beta_prior

        map_mcmc <- gMAP(cbind(mean, se) ~ 1 | study, data = hd_map,
                         family = gaussian, weights = n,
                         tau.dist = tau_dist, tau.prior = tau_p, beta.prior = beta_p)
        map <- automixfit(map_mcmc, type = "norm")
        if(!is.null(sigma)) sigma(map) <- sigma_ref else sigma_ref <- sigma(map)

        pred_dist  <- RBesT::preddist(map, n = n0, sigma = sigma_ref)
        pred_int   <- .pred_int_from(pred_dist, limits, alpha)
        map_robust <- robustify(map, weight = robust,
                                mean = summary(map)["mean"], sigma = sigma_ref)
        cap <- sum(hd$n)

        list(prior      = map,
             prior_rob  = map_robust,
             pred_dist  = pred_dist,
             pred_int   = pred_int,
             eff_n_map  = .ess_safe(map,        cap),
             eff_n      = .ess_safe(map_robust, cap),
             weights    = data.frame(w_prior = c(1 - robust, robust),
                                     w_post  = c(NA, NA)),
             tau_prior  = .tau_prior_out(tau_dist, tau_p),
             beta_prior = beta_p)
}

#-------------------------------------------------------------------------------
# Internal: multi-arm PI_norm(). One prediction interval per named arm. 'n0' is
# a single value (applied to all arms) or a named vector (one per arm).
#-------------------------------------------------------------------------------
.pi_norm_multiarm <- function(histdat, n0, robust, sigma,
                              tau_prior, tau_dist, beta_prior, limits, alpha, seed){

        if(!is.null(seed)) set.seed(seed)
        arms <- names(histdat)
        n0v  <- if(length(n0) == 1) stats::setNames(rep(n0, length(arms)), arms) else n0[arms]
        if(any(is.na(n0v)))
                stop("A named 'n0' must supply a value for every arm in 'histdat'.")

        b <- lapply(arms, function(a)
                    .pi_norm_arm(histdat[[a]], n0v[[a]], robust, sigma,
                                 tau_prior, tau_dist, beta_prior, limits, alpha))
        names(b) <- arms
        pick  <- function(nm) lapply(b, `[[`, nm)
        pickv <- function(nm) vapply(b, `[[`, numeric(1), nm)
        na_v  <- stats::setNames(rep(NA_real_, length(arms)), arms)

        out_list <- list("comp"        = NA,
                         "contr"       = NA,
                         "histdat"     = histdat,
                         "newdat"      = NA,
                         "distr"       = "normal",
                         "borrow_arms" = arms,
                         "prior_par"   = NA,
                         "weights"     = pick("weights"),
                         "eff_n_map"   = pickv("eff_n_map"),
                         "eff_n"       = pickv("eff_n"),
                         "eff_n_post"  = na_v,
                         "estimates"   = NA,
                         "contr_mat"   = NA,
                         "prior"       = pick("prior"),
                         "prior_rob"   = pick("prior_rob"),
                         "posterior"   = NA,
                         "pred_dist"   = pick("pred_dist"),
                         "pred_int"    = pick("pred_int"),
                         "shrinkage"   = NA,
                         "simpost"     = NA,
                         "tau_prior"   = pick("tau_prior"),
                         "beta_prior"  = pick("beta_prior"),
                         "alpha"       = alpha,
                         "seed"        = if(is.null(seed)) NA_integer_ else seed)

        structure(out_list, class = c("PI", "MultiMAP_norm", "MultiMAP"))
}

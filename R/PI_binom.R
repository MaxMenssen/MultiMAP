#' Prior and prediction interval for a binomial endpoint
#'
#' Builds robustified priors for the study arms of interest -- either a meta-analytic-predictive (MAP) prior
#' (\code{prior_est = "MAP"}, via \code{RBesT::gMAP})
#' or an empirical-Bayes beta prior estimated from the historical data
#' (\code{prior_est = "EB"}) -- and derives the prediction intervals for the
#' binomial proportions in the study arms of interest.
#'
#' @param histdat external (historical) data for the study-arms of interest; 
#' borrowing for one study arm (e.g. control group): submit a data.frame with columns
#' \code{events} and \code{non_events} (successes and failures). 
#' To borrow for more than one arm, pass a \emph{named list}
#'   containing different data frames: each element supplies the external
#'    data for a particular study-arm. Every arm without a list entry uses the 
#'    non-borrowing (flat-prior) posterior.
#' @param n0 planned sample size (e.g. number of patients) in the study arms of interest
#'   for which the prediction intervals are derived. 
#'   Either a single value (used for every arm) or a named vector with one value per arm.
#' @param prior_est prior for the concurrent control proportion:
#'   \code{"MAP"} (a meta-analytic-predictive prior via \code{RBesT::gMAP},
#'   \code{family = binomial}, fitted to a beta mixture) or \code{"EB"} (an
#'   empirical-Bayes beta prior whose parameters are estimated from the
#'   historical control data via \code{predint::pi_rho_est}, following
#'   Lui, 2000).
#' @param limits limits of the prediction interval: either \code{"two.sided"},
#'   \code{"lower"} or \code{"upper"}.
#' @param robust weight of the weakly-informative (robust) prior component added
#'   to robustify the prior.
#' @param tau_prior optional override for the \code{gMAP} between-study
#'   heterogeneity prior (\code{tau.prior} in \code{RBesT::gMAP}), on the logit
#'   scale. \code{NULL} (the default) uses the built-in heuristic (a
#'   \code{HalfNormal} scale of 1 for a mean proportion in [0.2, 0.8], otherwise
#'   \code{round(1 / sqrt(p_bar (1 - p_bar)), 1) / 2}); a larger scalar scale
#'   widens the MAP prior and its prediction interval when few historical studies
#'   leave tau weakly identified.
#' @param tau_dist distribution family of the \code{gMAP} between-study
#'   heterogeneity prior (\code{tau.dist} in \code{RBesT::gMAP}); one of
#'   \code{"HalfNormal"} (the default), \code{"TruncNormal"}, \code{"Uniform"},
#'   \code{"Gamma"}, \code{"InvGamma"}, \code{"LogNormal"}, \code{"TruncCauchy"},
#'   \code{"Exp"} or \code{"Fixed"}. The built-in default \code{tau_prior}
#'   heuristic is calibrated for the half-normal, so when \code{tau_dist} is set
#'   to any other family a matching \code{tau_prior} must also be supplied.
#' @param beta_prior optional override for the \code{gMAP} intercept prior
#'   (\code{beta.prior} in \code{RBesT::gMAP}), on the logit scale. \code{NULL}
#'   (the default) uses a scale of 2.
#' @param alpha 1 - coverage level of the prediction interval.
#' @param n_mcmc total number of Monte-Carlo samples to draw.
#' @param seed optional integer; if supplied, \code{set.seed(seed)} is called so
#'   the \code{gMAP} MCMC (for \code{prior_est = "MAP"}) is reproducible.
#'   \code{NULL} (default) leaves the RNG state untouched.
#'
#' @return A MultiMAP object: a named list with S3 class
#'   \code{c("PI", "MultiMAP_binom", "MultiMAP")} and the following slots:
#'   \describe{
#'     \item{\code{comp}}{\code{NA} (only \code{\link{SCI_binom}} fills this).}
#'     \item{\code{contr}}{\code{NA} (only \code{\link{SCI_binom}} fills this).}
#'     \item{\code{histdat}}{the historical control data used.}
#'     \item{\code{newdat}}{\code{NA} (no current trial groups in a
#'       prediction-interval fit).}
#'     \item{\code{distr}}{the endpoint family, \code{"betabinomial"}.}
#'     \item{\code{prior_par}}{\code{NA} (only \code{\link{SCI_binom}} fills this).}
#'     \item{\code{weights}}{data frame of the prior weights of the robustified
#'       mixture's components (posterior weights are \code{NA}).}
#'     \item{\code{eff_n_map}}{effective sample size of the (non-robustified)
#'       prior.}
#'     \item{\code{eff_n}}{effective sample size of the robustified prior.}
#'     \item{\code{eff_n_post}}{\code{NA} (no posterior is formed).}
#'     \item{\code{estimates}}{\code{NA} (only \code{\link{SCI_binom}} fills this).}
#'     \item{\code{contr_mat}}{\code{NA} (only \code{\link{SCI_binom}} fills this).}
#'     \item{\code{prior}}{the (non-robustified) prior (RBesT beta mixture).}
#'     \item{\code{prior_rob}}{the robustified prior (RBesT beta mixture).}
#'     \item{\code{posterior}}{\code{NA} (no posterior is formed).}
#'     \item{\code{pred_dist}}{the prior-predictive distribution for the number of
#'       events in a future concurrent control group, an RBesT mixture object
#'       built with \code{\link[RBesT]{preddist}}.}
#'     \item{\code{pred_int}}{the prediction interval for the concurrent control
#'       proportion (named \code{median}, \code{lower}, \code{upper}).}
#'     \item{\code{shrinkage}}{\code{NA} (only \code{\link{SCI_binom}} fills this).}
#'     \item{\code{simpost}}{\code{NA} (only \code{\link{SCI_binom}} fills this).}
#'     \item{\code{tau_prior}}{a list with the gMAP between-study heterogeneity
#'       prior distribution name (\code{dist}) and its parameters (\code{prior})
#'       actually used (\code{NA} for the EB prior).}
#'     \item{\code{beta_prior}}{the gMAP intercept prior actually used
#'       (\code{NA} for the EB prior).}
#'     \item{\code{alpha}}{the \code{alpha} used.}
#'     \item{\code{seed}}{the \code{seed} used (\code{NA} if none).}
#'   }
#'   Use \code{\link{summary.MultiMAP}}, \code{\link{get_output}} and
#'   \code{\link{report_binom}} on it.
#'
#' @seealso \code{\link{SCI_binom}}, \code{\link{report_binom}}
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
#' Lui, K.-J. (2000). A revisit of some simple point estimators in binomial
#' distributions. \emph{Biometrical Journal, 42}(6), 733--741.
#'
#' Zhao, Q., Zhu, Y., Xu, Z., Cheng, Z., Mei, J., Chen, X., & Wang, X. (2018).
#' Effect of ticagrelor plus aspirin, ticagrelor alone, or aspirin alone on
#' saphenous vein graft patency 1 year after coronary artery bypass grafting: a
#' randomized clinical trial (DACAB). \emph{JAMA, 319}(16), 1677--1686.
#' \doi{10.1001/jama.2018.3197}
#'
#' Chiarito, M., Sanz-Sanchez, J., Cannata, F., Cao, D., Sturla, M., Panico, C.,
#' Godino, C., Regazzoli, D., Reimers, B., De Caterina, R., Condorelli, G.,
#' Ferrante, G., & Stefanini, G. G. (2020). Monotherapy with a P2Y12 inhibitor or
#' aspirin for secondary prevention in patients with established atherosclerosis:
#' a systematic review and meta-analysis. \emph{Lancet, 395}(10235), 1487--1495.
#' \doi{10.1016/S0140-6736(20)30315-9}
#'
#' @examples
#' \donttest{
#' # Example data on pathological findings in rodent long-term carcinogenicity
#' # studies from EFSA 2025 as analysed in Menssen et al. 2026
#'
#' # HCD
#' histdat <- data.frame(events     = c(2, 2, 1, 5, 3, 5, 2, 2,
#'                                      4, 2, 3, 3, 3, 2, 3, 4),
#'                       non_events = c(48, 47, 49, 44, 51, 44, 54, 53,
#'                                      46, 46, 46, 86, 71, 45, 45, 46))
#'
#' # Compute the prediction interval for a future control group of 50 animals
#' m_pi <- PI_binom(histdat = histdat, n0 = 50, prior_est = "MAP",
#'                  seed = 81771)
#' summary(m_pi)
#' }
#'
#' \donttest{
#' # Antiplatelet therapy after coronary artery bypass grafting (DACAB trial,
#' # Zhao et al. 2018). External myocardial-infarction data are borrowed for the
#' # two monotherapy arms, P2Y12 (ticagrelor) and ASA (aspirin), from the
#' # Chiarito et al. 2020 meta-analysis. A named list of HCD gives one prediction
#' # interval per arm.
#'
#' # HCD: one data frame per borrowing arm (one row per historical study:
#' # CAPRIE, STAMI, AAASPS, CADET, ASCET, SOCRATES, TiCAB)
#' histdat <- list(
#'   P2Y12 = data.frame(events     = c(275,   8,   9,   1,  18,   25,  19),
#'                      non_events = c(9324, 726, 893,  93, 481, 6564, 912)),
#'   ASA   = data.frame(events     = c(333,  18,   8,   6,  18,   21,  30),
#'                      non_events = c(9253, 718, 899,  84, 484, 6589, 898)))
#'
#'# Named vector for different sample size
#'n_0 <- c("P2Y12"=166, "ASA"=167)
#'
#'# Prediction intervals for the two arms of interest
#'m_pi <- PI_binom(histdat = histdat, n0 = n_0, prior_est = "MAP", seed = 2024)
#'
#' summary(m_pi)
#' }
#' @export
PI_binom <- function(histdat,
                     n0,
                     prior_est = "MAP",
                     limits = "two.sided",
                     robust = 0.2,
                     tau_prior = NULL,
                     tau_dist = "HalfNormal",
                     beta_prior = NULL,
                     alpha = 0.05,
                     n_mcmc = 50000,
                     seed = NULL){

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
                stop("'histdat' must be a data frame with at least 2 columns (events, non-events).")
        if(!is.numeric(n0) || !all(n0 > 0))
                stop("'n0' must be positive (a single value or one per borrowed arm).")
        if(!multi_arm && length(n0) != 1)
                stop("'n0' must be a single positive number.")
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
                return(.pi_binom_multiarm(histdat, n0, prior_est, robust,
                                          tau_prior, tau_dist, beta_prior, limits, alpha, seed))
        }

        # Effective gMAP hyperpriors actually used (filled in the MAP block; NA
        # for the EB path). Stored in the output.
        tau_p  <- NA
        beta_p <- NA

        #=======================================================================
        # Prior specification -- same as SCI_binom(): MAP (gMAP) or EB (emp. Bayes)
        #=======================================================================

        if(prior_est=="MAP"){

                histdat_map <- cbind(data.frame(study=factor(1:nrow(histdat))),
                                     histdat)

                # Default tau.prior (HalfNormal scale on the logit scale) depends
                # on the mean proportion p_bar.
                p_bar <- mean(histdat_map[,2] / (histdat_map[,2] + histdat_map[,3]))

                if(0.2 <= p_bar & p_bar <= 0.8){
                        tau_default <- 1
                }
                else{
                        if(all(histdat_map[,2]==0)){
                                y_vec <- histdat_map[,2]
                                y_vec[1] <- 1
                                n_vec <- histdat_map[,2] + histdat_map[,3]
                                n_vec[1] <- n_vec[1] -1
                                p_bar <- mean(y_vec / (n_vec))
                        }
                        tau_default <- round(1 / sqrt(p_bar * (1 - p_bar)), 1) / 2
                }

                # tau_prior / beta_prior override the defaults when supplied.
                tau_p  <- if(is.null(tau_prior))  tau_default else tau_prior
                beta_p <- if(is.null(beta_prior)) 2           else beta_prior

                map_mcmc <- gMAP(cbind(histdat_map[,2], histdat_map[,3]) ~ 1 | study,
                                 data = histdat_map,
                                 tau.dist = tau_dist,
                                 tau.prior = tau_p,
                                 beta.prior = beta_p,
                                 family = binomial)

                # Derive beta mixture
                map <- automixfit(map_mcmc)
        }

        #-----------------------------------------------------------------------

        if(prior_est=="EB"){

                # No hist. observations (set a and b to predefined values)
                if(all(histdat[,1]==0)){
                        a_est <- 1
                        b_est <- sum(histdat[,2])-1
                }

                else{
                        # Get Estimates from HCD
                        estimates <- predint::pi_rho_est(histdat)
                        names(estimates) <- c("pi", "rho")
                        pi_HCD <- unname(estimates[1])
                        rho_HCD <- max(0.00001, unname(estimates[2]))

                        if(rho_HCD==0.00001){
                                warning("Data is underdispersed. rho_HCD was set to 0.00001")
                        }

                        # Estimates for a and b
                        a_b_HCD <- (1-rho_HCD)/rho_HCD
                        a_est <- pi_HCD*a_b_HCD
                        b_est <- a_b_HCD-a_est

                        # adjust a and b if data is not overdispersed
                        if(a_est > sum(histdat[,1]) & b_est > sum(histdat[,1])){
                                a_est <- sum(histdat[,1])
                                b_est <- sum(histdat[,2])
                        }
                }

                map <- mixbeta(c(1, a_est, b_est))
        }

        #=======================================================================
        # Robustify + prediction interval
        #=======================================================================

        # Add weakly-informative (robust) component (centred at 1/2).
        map_robust <- robustify(map, weight = robust, mean = 1 / 2)

        #-----------------------------------------------------------------------
        # Prediction interval for the concurrent control proportion.
        # preddist() on the (non-robustified) prior gives the beta-binomial
        # predictive for the number of events in n0 future control patients;
        # dividing the count quantiles by n0 puts the interval on the
        # proportion scale.

        pred_dist <- RBesT::preddist(map, n = n0)

        if(limits == "two.sided"){
                pred_int <- c(median = unname(RBesT::qmix(pred_dist, 0.5))       / n0,
                              lower  = unname(RBesT::qmix(pred_dist, alpha/2))    / n0,
                              upper  = unname(RBesT::qmix(pred_dist, 1 - alpha/2))/ n0)
        }

        if(limits == "lower"){
                pred_int <- c(median = unname(RBesT::qmix(pred_dist, 0.5))  / n0,
                              lower  = unname(RBesT::qmix(pred_dist, alpha)) / n0,
                              upper  = NA)
        }

        if(limits == "upper"){
                pred_int <- c(median = unname(RBesT::qmix(pred_dist, 0.5))       / n0,
                              lower  = NA,
                              upper  = unname(RBesT::qmix(pred_dist, 1 - alpha)) / n0)
        }

        #-----------------------------------------------------------------------
        # Effective sample size of the (non-robustified) prior

        eff_n_map <- try(round(ess(map, method = "elir")), silent = TRUE)

        if(inherits(eff_n_map, "try-error")){

                eff_n_map <- round(ess(map, method = "moment"))

                warning("At least one parameter of the beta mixtures is less than 1.\n Estimation of ESS via elir-method does not work. ESS was estimated based on method of moments")
        }

        eff_n_map <- min(max(0, eff_n_map),
                         sum(c(histdat[,1] + histdat[,2])))

        #-----------------------------------------------------------------------
        # Effective sample size of the robustified prior

        eff_sample_size <- try(round(ess(map_robust, method = "elir")), silent = TRUE)

        if(inherits(eff_sample_size, "try-error")){

                eff_sample_size <- round(ess(map_robust, method = "moment"))

                warning("At least one parameter of the beta mixtures is less than 1.\n Estimation of ESS via elir-method does not work. ESS was estimated based on method of moments")
        }

        eff_sample_size <- min(max(0, eff_sample_size),
                               sum(c(histdat[,1] + histdat[,2])))

        #-----------------------------------------------------------------------
        weights <- data.frame(w_prior = c(1 - robust, robust),
                              w_post  = c(NA, NA))

        #-----------------------------------------------------------------------
        # Output

        out_list <- list("comp"       = NA,
                         "contr"      = NA,
                         "histdat"    = histdat,
                         "newdat"     = NA,
                         "distr"      = "betabinomial",
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
                            class = c("PI", "MultiMAP_binom", "MultiMAP"))

        return(out_s3)
}


#-------------------------------------------------------------------------------
# Internal: MAP/EB prior + prediction interval for ONE arm (binomial, no
# posterior update). Reuses .map_prior_binom() (defined in SCI_binom.R).
#-------------------------------------------------------------------------------
.pi_binom_arm <- function(hd, n0, prior_est, robust, tau_prior, tau_dist, beta_prior, limits, alpha){

        pr  <- .map_prior_binom(hd, prior_est, tau_prior, tau_dist, beta_prior)
        map <- pr$map
        map_robust <- robustify(map, weight = robust, mean = 1 / 2)
        pred_dist  <- RBesT::preddist(map, n = n0)
        pred_int   <- .pred_int_from(pred_dist, limits, alpha, scale = n0)
        cap        <- sum(hd[, 1] + hd[, 2])

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
# Internal: multi-arm PI_binom(). One prediction interval per named arm.
#-------------------------------------------------------------------------------
.pi_binom_multiarm <- function(histdat, n0, prior_est, robust,
                               tau_prior, tau_dist, beta_prior, limits, alpha, seed){

        if(!is.null(seed)) set.seed(seed)
        arms <- names(histdat)
        n0v  <- if(length(n0) == 1) stats::setNames(rep(n0, length(arms)), arms) else n0[arms]
        if(any(is.na(n0v)))
                stop("A named 'n0' must supply a value for every arm in 'histdat'.")

        b <- lapply(arms, function(a)
                    .pi_binom_arm(histdat[[a]], n0v[[a]], prior_est, robust,
                                  tau_prior, tau_dist, beta_prior, limits, alpha))
        names(b) <- arms
        pick  <- function(nm) lapply(b, `[[`, nm)
        pickv <- function(nm) vapply(b, `[[`, numeric(1), nm)
        na_v  <- stats::setNames(rep(NA_real_, length(arms)), arms)

        out_list <- list("comp"        = NA, "contr" = NA, "histdat" = histdat,
                         "newdat"      = NA, "distr" = "betabinomial",
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

        structure(out_list, class = c("PI", "MultiMAP_binom", "MultiMAP"))
}

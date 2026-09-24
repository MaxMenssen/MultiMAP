#' Bayesian dynamic borrowing with simultaneous credible intervals (binomial endpoint)
#'
#' Bayesian dynamic borrowing for a binomial endpoint following Menssen et al. 2026. 
#' The pipeline: (1) build informative priors for the study arms of interest
#' based on external data (such as historical controls) -- either meta-analytic-predictive (MAP)
#' or empirical-Bayes beta priors (2) robustify priors by mixing in a
#' non-informative beta(1,1) component; (3) update to posteriors with the current
#' control data; (4) draw samples from the joint posterior; 
#' (5) derive SIMULTANEOUS credible sets for the comparisons of interest 
#' based on the joint posterior 
#'
#' @param histdat external (historical) data for the study-arms of interest; 
#' borrowing for one study arm (e.g. control group): submit a data.frame with columns
#' \code{events} and \code{non_events} (successes and failures). 
#' To borrow for more than one arm, pass a \emph{named list}
#'   containing different data frames: each element supplies the external
#'    data for a particular study-arm. 
#' @param newdat a data.frame containing the current trial; columns \code{group}, \code{events},
#'   \code{non_events}. In case for borrowing for the control group: Row 1 is the concurrent 
#'   control, the remaining rows are treatment groups.
#' @param prior_est \code{"MAP"} (a meta-analytic-predictive prior via \code{RBesT::gMAP}) 
#'   or \code{"EB"} (empirical-Bayes beta prior).
#' @param type type of comparison: either
#'   \code{"Dunnett"} (many-to-one comparisons) or
#'   \code{"Tukey"} (all pairwise comparisons).
#' @param contr  \code{"OR"} provides odds ratios; \code{"RR"} provides ratios of proportions (risk-ratios);
#' \code{"pi_diff"} provides differences of proportions.
#' @param limits limits of the credible intervals: either \code{"two.sided"}, \code{"lower"} or \code{"upper"}.
#' @param limits_pi limits of the prediction interval for the concurrent
#'   control proportion: either \code{"two.sided"}, \code{"lower"} or
#'   \code{"upper"}. Defaults to \code{limits} so that, unless set, the
#'   prediction interval follows the contrasts' one- or two-sidedness.
#' @param base index of the control group for the contrast matrix.
#' @param robust weight of the weakly-informative (robust) prior component.
#' @param tau_prior optional override for the \code{gMAP} between-study
#'   heterogeneity prior (\code{tau.prior} in \code{RBesT::gMAP}), on the logit
#'   scale. If \code{NULL}, a \code{HalfNormal} scale of 1 for a mean proportion 
#'   in [0.2, 0.8] is used, otherwise the \code{HalfNormal} scale is set to 
#'   \code{(1 / sqrt(p_bar (1 - p_bar)) / 2}). Supply any \code{cbind(0, scale)} to override.
#' @param tau_dist distribution family of the \code{gMAP} between-study
#'   heterogeneity prior (\code{tau.dist} in \code{RBesT::gMAP}); one of
#'   \code{"HalfNormal"} (the default), \code{"TruncNormal"}, \code{"Uniform"},
#'   \code{"Gamma"}, \code{"InvGamma"}, \code{"LogNormal"}, \code{"TruncCauchy"},
#'   \code{"Exp"} or \code{"Fixed"}. The built-in default \code{tau_prior}
#'   heuristic is calibrated for the half-normal, so when \code{tau_dist} is set
#'   to any other family a matching \code{tau_prior} must also be supplied.
#' @param beta_prior optional override for the \code{gMAP} intercept prior
#'   (\code{beta.prior} in \code{RBesT::gMAP}), on the logit scale. If \code{NULL},
#'   a normal(0, 2) is used.
#' @param alpha 1 - simultaneous confidence level.
#' @param n_mcmc total number of Monte-Carlo samples to draw.
#' @param seed optional integer; if supplied, \code{set.seed(seed)} is called so
#'   the \code{gMAP} MCMC (for \code{prior_est = "MAP"}) and the Monte-Carlo
#'   sampling are reproducible. \code{NULL} (default) leaves the RNG state
#'   untouched.
#'
#' @return A MultiMAP object: a named list with S3 class
#'   \code{c("SCI", "MultiMAP_binom", "MultiMAP")} and the following slots:
#'   \describe{
#'     \item{\code{comp}}{data frame of the simultaneous credible intervals: one
#'       row per contrast with columns \code{comp} (label), \code{estimate},
#'       \code{lower.ci}, \code{upper.ci}.}
#'     \item{\code{contr}}{the contrast type used (\code{"RR"}, \code{"OR"}
#'       or \code{"pi_diff"}).}
#'     \item{\code{histdat}}{the historical control data used.}
#'     \item{\code{newdat}}{the current trial data (concurrent control in row 1),
#'       with the observed per-group proportion and posterior beta parameters
#'       added.}
#'     \item{\code{distr}}{the endpoint family, \code{"betabinomial"}.}
#'     \item{\code{prior_par}}{the beta parameters (\code{a}, \code{b}) of the
#'       robustified prior mixture's components.}
#'     \item{\code{weights}}{data frame of the prior and posterior weights of the
#'       robustified mixture's components.}
#'     \item{\code{eff_n_map}}{effective sample size of the (non-robustified)
#'       prior.}
#'     \item{\code{eff_n}}{effective sample size of the robustified prior.}
#'     \item{\code{eff_n_post}}{effective sample size of the concurrent-control
#'       posterior.}
#'     \item{\code{estimates}}{the internal sampling table (mixture weights, beta
#'       parameters and per-component draw counts).}
#'     \item{\code{contr_mat}}{the contrast matrix used.}
#'     \item{\code{prior}}{the (non-robustified) prior (RBesT beta mixture).}
#'     \item{\code{prior_rob}}{the robustified prior (RBesT beta mixture).}
#'     \item{\code{posterior}}{the concurrent-control posterior (RBesT beta
#'       mixture).}
#'     \item{\code{pred_dist}}{the prior-predictive distribution for the number of
#'       events in a future concurrent control group.}
#'     \item{\code{pred_int}}{the prediction interval for the concurrent control
#'       proportion (named \code{median}, \code{lower}, \code{upper}).}
#'     \item{\code{shrinkage}}{named numeric \code{c(abs, rel)}: the absolute
#'       shrinkage of the concurrent control proportion towards the borrowed
#'       prior and its relative version (absolute shrinkage divided by the raw
#'       concurrent control proportion).}
#'     \item{\code{simpost}}{matrix of posterior Monte-Carlo samples (control plus
#'       one column per treatment group).}
#'     \item{\code{tau_prior}}{a list with the gMAP between-study heterogeneity
#'       prior distribution name (\code{dist}) and its parameters (\code{prior})
#'       actually used (\code{NA} for the EB prior).}
#'     \item{\code{beta_prior}}{the gMAP intercept prior actually used
#'       (\code{NA} for the EB prior).}
#'     \item{\code{alpha}}{the \code{alpha} used.}
#'     \item{\code{seed}}{the \code{seed} used (\code{NA} if none).}
#'   }
#'
#' @seealso \code{\link{SCI_norm}}
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
#' Besag, J., Green, P., Higdon, D., & Mengersen, K. (1995). Bayesian computation
#' and stochastic systems. \emph{Statistical Science, 10}(1), 3--66.
#' \url{http://www.jstor.org/stable/2246224}
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
#' ### One arm: Borrowing from historical control data
#' 
#' # Example data on pathological findings in rodent long-term carcinogenicity
#' # studies from EFSA 2025 as analysed in Menssen et al. 2026
#' 
#' # HCD
#' histdat <- data.frame(events     = c(2, 2, 1, 5, 3, 5, 2, 2, 
#'                                      4, 2, 3, 3, 3, 2, 3, 4),
#'                       non_events = c(48, 47, 49, 44, 51, 44, 54, 53,
#'                                      46, 46, 46, 86, 71, 45, 45, 46))
#'
#' # Current trial
#' newdat  <- data.frame(group      = c("0", "10", "50", "1500", "10000"),
#'                       events     = c(3, 4, 5, 7, 8),
#'                       non_events = c(44, 36, 39, 42, 41))
#'
#' # Calculate simultaneous credible intervals
#' res <- SCI_binom(histdat, newdat, prior_est = "MAP",
#'                  contr = "RR", seed = 81771)
#' summary(res)
#' }
#'
#' \donttest{
#' #----------------------------------------------------------------------------
#' ### Multi-arm borrowing
#' 
#' # Antiplatelet therapy after coronary artery bypass grafting (DACAB trial,
#' # Zhao et al. 2018). External myocardial-infarction data are borrowed for the
#' # two monotherapy arms, P2Y12 (ticagrelor) and ASA (aspirin), from the
#' # Chiarito et al. 2020 meta-analysis; the DAPT combination arm has no external
#' # data and is the reference. A named list of HCD triggers multi-arm borrowing.
#'
#' # HCD: one data frame per borrowing arm (one row per historical study:
#' # CAPRIE, STAMI, AAASPS, CADET, ASCET, SOCRATES, TiCAB)
#' histdat <- list(
#'   P2Y12 = data.frame(events     = c(275,   8,   9,   1,  18,   25,  19),
#'                      non_events = c(9324, 726, 893,  93, 481, 6564, 912)),
#'   ASA   = data.frame(events     = c(333,  18,   8,   6,  18,   21,  30),
#'                      non_events = c(9253, 718, 899,  84, 484, 6589, 898)))
#'
#' # Current trial: DAPT is the reference (row 1), P2Y12 and ASA are borrowed
#' newdat  <- data.frame(group      = c("DAPT", "P2Y12", "ASA"),
#'                       events     = c(2, 2, 3),
#'                       non_events = c(166, 164, 163))
#'
#' # Simultaneous credible intervals for the odds ratios vs DAPT
#' res <- SCI_binom(histdat, newdat, prior_est = "MAP",
#'                  type = "Dunnett", contr = "OR", base = 1, seed = 2024)
#' summary(res)
#' }
#' @export
SCI_binom <- function(histdat,
                      newdat,
                      prior_est = "MAP",
                      type = "Dunnett",
                      contr = "OR",
                      limits = "two.sided",
                      limits_pi = limits,
                      base = 1,
                      robust = 0.2,
                      tau_prior = NULL,
                      tau_dist = "HalfNormal",
                      beta_prior = NULL,
                      alpha = 0.05,
                      n_mcmc = 50000,
                      seed = NULL){

        # Optional reproducibility: seeds gMAP()'s MCMC (MAP path) and the
        # Monte-Carlo sampling below.
        if(!is.null(seed)) set.seed(seed)

        # User-facing one-sidedness vocabulary: two.sided / lower / upper.
        # 'limits' drives the contrasts, 'limits_pi' the prediction
        # interval (independent; defaults to 'limits').
        limits    <- match.arg(limits,    c("two.sided", "lower", "upper"))
        limits_pi <- match.arg(limits_pi, c("two.sided", "lower", "upper"))

        # Contrast type: ratio of proportions / odds ratio / difference of props.
        contr <- match.arg(contr, c("RR", "OR", "pi_diff"))

        # Multiple-comparison type: many-to-one (Dunnett) or all-pairwise (Tukey).
        type <- match.arg(type, c("Dunnett", "Tukey"))

        # gMAP heterogeneity-prior distribution (tau.dist); non-half-normal
        # families need an explicit tau_prior.
        tau_dist <- match.arg(tau_dist, .tau_dist_choices)
        .check_tau_dist(tau_dist, tau_prior)

        # Multi-arm borrowing when 'histdat' is a named list (one historical data
        # set per arm); a plain data frame keeps the single control-only behaviour.
        multi_arm <- .is_multi_arm(histdat)

        #-----------------------------------------------------------------------
        # Input checks

        if(!multi_arm && (!is.data.frame(histdat) || ncol(histdat) < 2))
                stop("'histdat' must be a data frame with at least 2 columns (events, non-events).")
        if(!is.data.frame(newdat) || ncol(newdat) < 3)
                stop("'newdat' must be a data frame with at least 3 columns (group, events, non-events).")
        if(!prior_est %in% c("MAP", "EB"))
                stop("'prior_est' must be either 'MAP' or 'EB'.")
        if(!is.numeric(robust) || length(robust) != 1 || robust < 0 || robust > 1)
                stop("'robust' must be a single number in [0, 1].")
        if(!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1)
                stop("'alpha' must be a single number in (0, 1).")
        if(!is.numeric(base) || length(base) != 1 || base < 1 || base > nrow(newdat))
                stop("'base' must be a single group index in 1:nrow(newdat).")

        #-----------------------------------------------------------------------
        # Multi-arm borrowing: each named arm gets its own prior + posterior; the
        # remaining arms use the non-borrowing Beta(1 + ev, 1 + non-ev) posterior.
        if(multi_arm){
                .check_multi_histdat(histdat, newdat[[1]])
                return(.sci_binom_multiarm(histdat, newdat, prior_est, contr, type, base,
                                           robust, tau_prior, tau_dist, beta_prior,
                                           limits, limits_pi, alpha, n_mcmc, seed))
        }

        # Effective gMAP hyperpriors actually used (filled in the MAP block; NA
        # for the EB path, which uses no gMAP). Stored in the output.
        tau_p  <- NA
        beta_p <- NA

        #=======================================================================
        # Prior specification -- unchanged: MAP (gMAP) or EB (empirical Bayes)
        #=======================================================================

        if(prior_est=="MAP"){

                histdat_map <- cbind(data.frame(study=factor(1:nrow(histdat))),
                                     histdat)

                # Default tau.prior (HalfNormal scale on the logit scale) depends
                # on the mean proportion p_bar: 1 for a moderate p_bar in
                # [0.2, 0.8], otherwise the n_infinity heuristic for extreme rates.
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

                #---------------------------------------------------------------
                # Derive beta mixture
                map <- automixfit(map_mcmc)
        }


        #-----------------------------------------------------------------------

        if(prior_est=="EB"){

                # No hist. observations (set a and b to predefined values)
                if(all(histdat[,1]==0)){
                        a_est <- 1
                        b_est <- sum(histdat[,2])-1

                        estimates_out <- NULL
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

                        estimates_out <- c("pi"=pi_HCD, "rho"=rho_HCD)
                }

                map <- mixbeta(c(1, a_est, b_est))
        }

        #=======================================================================
        # Robustify, update, sample
        #=======================================================================

        # Add uninformative prior
        map_robust <- robustify(map, weight = robust, mean = 1 / 2)

        map_robust_df <- data.frame(cbind(t(map_robust)))
        colnames(map_robust_df) <- c("w", "a", "b")

        #-----------------------------------------------------------------------
        # Prediction interval for the concurrent control proportion.
        # preddist() on the (non-robustified) prior gives the beta-binomial
        # predictive for the number of events in n0 future control patients;
        # dividing the count quantiles by n0 puts the interval on the
        # proportion scale.

        n0 <- newdat[1, 2] + newdat[1, 3]

        pred_dist <- RBesT::preddist(map, n = n0)

        if(limits_pi == "two.sided"){
                pred_int <- c(median = unname(RBesT::qmix(pred_dist, 0.5))       / n0,
                              lower  = unname(RBesT::qmix(pred_dist, alpha/2))    / n0,
                              upper  = unname(RBesT::qmix(pred_dist, 1 - alpha/2))/ n0)
        }

        if(limits_pi == "lower"){
                pred_int <- c(median = unname(RBesT::qmix(pred_dist, 0.5))  / n0,
                              lower  = unname(RBesT::qmix(pred_dist, alpha)) / n0,
                              upper  = NA)
        }

        if(limits_pi == "upper"){
                pred_int <- c(median = unname(RBesT::qmix(pred_dist, 0.5))       / n0,
                              lower  = NA,
                              upper  = unname(RBesT::qmix(pred_dist, 1 - alpha)) / n0)
        }

        # Translate the user-facing 'limits' (two.sided / lower / upper) to
        # the BSagri::SCSnp() vocabulary (two.sided / greater / less): a lower
        # simultaneous bound corresponds to 'greater', an upper bound to 'less'.
        alt_scs <- switch(limits,
                          "lower" = "greater",
                          "upper" = "less",
                          "two.sided")

        #-----------------------------------------------------------------------
        # Posterior for the concurrent control group

        posterior <- postmix(map_robust, n=newdat[1,2] + newdat[1,3], r=newdat[1,2])

        post_df <- data.frame(cbind(t(posterior)))
        colnames(post_df) <- c("w_post", "post_a", "post_b")
        post_df <- post_df[, c("post_a", "post_b", "w_post")]

        simdat <- data.frame(cbind(map_robust_df, post_df))
        simdat$n       <- n_mcmc                          # total no. of samples to draw
        simdat$n_prior <- round(simdat$w      * simdat$n) # samples from each prior comp.
        simdat$n_post  <- round(simdat$w_post * simdat$n) # samples from each posterior comp.

        #-----------------------------------------------------------------------
        # Sample from the beta posterior mixture (control group)

        simpost_ccg_fun <- function(x){
                rbeta(n      = x$n_post,
                      shape1 = x$post_a,
                      shape2 = x$post_b)
        }

        simpost_ccg <- unlist(lapply(X   = split(simdat, seq(nrow(simdat))),
                                     FUN = simpost_ccg_fun))

        #-----------------------------------------------------------------------
        # Posterior for the treatment groups
        # Beta(1 + events, 1 + non-events) with a uniform Beta(1, 1) prior.

        newdat$posta <- 1 + newdat[,2]
        newdat$postb <- 1 + newdat[,3]

        simpost_trt <- apply(X      = newdat[-1, c("posta", "postb")],
                             MARGIN = 1,
                             FUN    = function(x){rbeta(n      = length(simpost_ccg),
                                                        shape1 = x["posta"],
                                                        shape2 = x["postb"])})

        #-----------------------------------------------------------------------
        # Matrix with samples from posteriors

        simpost <- cbind(simpost_ccg,
                         simpost_trt)

        #-----------------------------------------------------------------------
        # Output object for weights

        weights <- data.frame(w_prior = c(1 - robust, robust),
                              w_post  = c(1 - simdat$w_post[nrow(simdat)],
                                          simdat$w_post[nrow(simdat)]))

        #-----------------------------------------------------------------------
        # Prior parameters + estimates output

        a_est <- simdat$a
        b_est <- simdat$b

        estimates_out <- simdat

        # Observed per-group proportion (concurrent control in row 1)
        newdat$mean <- newdat[, 2] / (newdat[, 2] + newdat[, 3])

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
        # Effective sample size of the concurrent-control posterior

        eff_n_post <- try(round(ess(posterior, method = "elir")), silent = TRUE)

        if(inherits(eff_n_post, "try-error")){

                eff_n_post <- round(ess(posterior, method = "moment"))

                warning("At least one parameter of the beta mixtures is less than 1.\n Estimation of ESS via elir-method does not work. ESS was estimated based on method of moments")
        }

        eff_n_post <- min(max(0, eff_n_post),
                          sum(c(histdat[,1] + histdat[,2])))

        #-----------------------------------------------------------------------
        # Shrinkage

        postmean <- summary(posterior)["mean"]
        ccg_mean <- newdat[1,2] / (newdat[1,2] + newdat[1,3])

        shrink_abs <- unname(abs(ccg_mean - postmean))
        shrink_rel <- if(ccg_mean == 0) NA_real_ else shrink_abs / abs(ccg_mean)
        shrinkage  <- c(abs = shrink_abs, rel = shrink_rel)

        #-----------------------------------------------------------------------
        ## ratio of proportions (pi_1 / pi_2)

        if(contr=="RR"){
                # Get contrast matrix
                nni <- newdat[,2] + newdat[,3]
                names(nni) <- newdat[,1]

                cmatr <- mratios::contrMatRatio(n=nni,
                                                type=type,
                                                base=base)

                #---------------------------------------------------------------
                # Compute ratios
                ratio2control <- BSagri::CCRatio.default(x=simpost, cmat=cmatr)
                dpostratio <- as.data.frame(ratio2control$chains)

                #---------------------------------------------------------------
                # Simultaneous confidence sets from empirical joint distribution
                SCIR <- BSagri::SCSnp(x=ratio2control, conf.level=1-alpha,
                                      alternative=alt_scs)

                dscir <- data.frame(
                        row.names = seq_along(SCIR$estimate),
                        comp      = names(SCIR$estimate),
                        estimate  = unname(signif(SCIR$estimate, 3)),
                        lower.ci  = unname(signif(SCIR$conf.int[,1], 3)),
                        upper.ci  = unname(signif(SCIR$conf.int[,2], 3)))
        }

        #-----------------------------------------------------------------------
        ## odds ratio

        if(contr=="OR"){
                # Get contrast matrix
                nni <- newdat[,2] + newdat[,3]
                names(nni) <- newdat[,1]

                cmatr <- multcomp::contrMat(n=nni,
                                            type=type,
                                            base=base)

                #---------------------------------------------------------------

                simpost <- log(simpost/(1-simpost))

                # Compute differences on the logit scale
                diff2control <- BSagri::CCDiff.default(x=simpost, cmat=cmatr)
                dpostdiff <- as.data.frame(diff2control$chains)

                #---------------------------------------------------------------
                # Simultaneous confidence sets from empirical joint distribution
                SCIR <- BSagri::SCSnp(x=diff2control, conf.level=1-alpha,
                                      alternative=alt_scs)

                dscir <- data.frame(
                        row.names = seq_along(SCIR$estimate),
                        comp      = names(SCIR$estimate),
                        estimate  = unname(signif(exp(SCIR$estimate), 3)),
                        lower.ci  = unname(signif(exp(SCIR$conf.int[,1]), 3)),
                        upper.ci  = unname(signif(exp(SCIR$conf.int[,2]), 3)))
        }

        #-----------------------------------------------------------------------
        ## difference of proportions (pi_1 - pi_2)

        if(contr=="pi_diff"){
                # Get contrast matrix
                nni <- newdat[,2] + newdat[,3]
                names(nni) <- newdat[,1]

                cmatr <- multcomp::contrMat(n=nni,
                                            type=type,
                                            base=base)

                #---------------------------------------------------------------
                # Compute differences
                diff2control <- BSagri::CCDiff.default(x=simpost, cmat=cmatr)
                dpostdiff <- as.data.frame(diff2control$chains)

                #---------------------------------------------------------------
                # Simultaneous confidence sets from empirical joint distribution
                SCIR <- BSagri::SCSnp(x=diff2control, conf.level=1-alpha,
                                      alternative=alt_scs)

                dscir <- data.frame(
                        row.names = seq_along(SCIR$estimate),
                        comp      = names(SCIR$estimate),
                        estimate  = unname(signif(SCIR$estimate, 3)),
                        lower.ci  = unname(signif(SCIR$conf.int[,1], 3)),
                        upper.ci  = unname(signif(SCIR$conf.int[,2], 3)))
        }

        #-----------------------------------------------------------------------
        # Output

        out_list <- list("comp"       = dscir,
                         "contr"      = contr,
                         "histdat"    = histdat,
                         "newdat"     = newdat,
                         "distr"      = "betabinomial",
                         "prior_par"  = c("a" = a_est, "b" = b_est),
                         "weights"    = weights,
                         "eff_n_map"  = eff_n_map,
                         "eff_n"      = eff_sample_size,
                         "eff_n_post" = eff_n_post,
                         "estimates"  = estimates_out,
                         "contr_mat"  = cmatr,
                         "prior"      = map,
                         "prior_rob"  = map_robust,
                         "posterior"  = posterior,
                         "pred_dist"  = pred_dist,
                         "pred_int"   = pred_int,
                         "shrinkage"  = shrinkage,
                         "simpost"    = simpost,
                         "tau_prior"  = .tau_prior_out(tau_dist, tau_p),
                         "beta_prior" = beta_p,
                         "alpha"      = alpha,
                         "seed"       = if(is.null(seed)) NA_integer_ else seed)

        out_s3 <- structure(out_list,
                            class = c("SCI", "MultiMAP_binom", "MultiMAP"))

        return(out_s3)
}


#-------------------------------------------------------------------------------
# Internal: build the informative prior (RBesT beta mixture) for one arm from its
# historical data (events, non-events), via MAP (gMAP) or EB (empirical Bayes).
# Shared by SCI_binom() and PI_binom(). Returns list(map, tau_p, beta_p, estimates).
#-------------------------------------------------------------------------------
.map_prior_binom <- function(hd, prior_est, tau_prior, tau_dist, beta_prior){

        tau_p <- NA; beta_p <- NA; estimates_out <- NA

        if(prior_est == "MAP"){
                hd_map <- cbind(data.frame(study = factor(seq_len(nrow(hd)))), hd)
                p_bar  <- mean(hd_map[, 2] / (hd_map[, 2] + hd_map[, 3]))
                if(0.2 <= p_bar & p_bar <= 0.8){
                        tau_default <- 1
                } else {
                        if(all(hd_map[, 2] == 0)){
                                y_vec <- hd_map[, 2]; y_vec[1] <- 1
                                n_vec <- hd_map[, 2] + hd_map[, 3]; n_vec[1] <- n_vec[1] - 1
                                p_bar <- mean(y_vec / n_vec)
                        }
                        tau_default <- round(1 / sqrt(p_bar * (1 - p_bar)), 1) / 2
                }
                tau_p  <- if(is.null(tau_prior))  tau_default else tau_prior
                beta_p <- if(is.null(beta_prior)) 2           else beta_prior

                map_mcmc <- gMAP(cbind(hd_map[, 2], hd_map[, 3]) ~ 1 | study,
                                 data = hd_map, tau.dist = tau_dist,
                                 tau.prior = tau_p, beta.prior = beta_p, family = binomial)
                map <- automixfit(map_mcmc)
        } else {
                # Empirical Bayes (Lui, 2000)
                if(all(hd[, 1] == 0)){
                        a_est <- 1; b_est <- sum(hd[, 2]) - 1
                } else {
                        est     <- predint::pi_rho_est(hd)
                        pi_HCD  <- unname(est[1])
                        rho_HCD <- max(0.00001, unname(est[2]))
                        if(rho_HCD == 0.00001)
                                warning("Data is underdispersed. rho_HCD was set to 0.00001")
                        a_b_HCD <- (1 - rho_HCD) / rho_HCD
                        a_est   <- pi_HCD * a_b_HCD
                        b_est   <- a_b_HCD - a_est
                        if(a_est > sum(hd[, 1]) & b_est > sum(hd[, 1])){
                                a_est <- sum(hd[, 1]); b_est <- sum(hd[, 2])
                        }
                        estimates_out <- c("pi" = pi_HCD, "rho" = rho_HCD)
                }
                map <- mixbeta(c(1, a_est, b_est))
        }
        list(map = map, tau_p = tau_p, beta_p = beta_p, estimates = estimates_out)
}

#-------------------------------------------------------------------------------
# Internal: borrowing pipeline for ONE arm (binomial endpoint). (ev_cur,
# nonev_cur) is the arm's current-trial count. Returns the per-arm bundle with
# n_mcmc posterior draws (proportions) in $samples.
#-------------------------------------------------------------------------------
.borrow_binom_arm <- function(hd, ev_cur, nonev_cur, prior_est, robust,
                              tau_prior, tau_dist, beta_prior, limits_pi, alpha, n_mcmc){

        pr  <- .map_prior_binom(hd, prior_est, tau_prior, tau_dist, beta_prior)
        map <- pr$map
        map_robust <- robustify(map, weight = robust, mean = 1 / 2)

        n0        <- ev_cur + nonev_cur
        pred_dist <- RBesT::preddist(map, n = n0)
        pred_int  <- .pred_int_from(pred_dist, limits_pi, alpha, scale = n0)

        posterior <- postmix(map_robust, n = n0, r = ev_cur)
        samples   <- RBesT::rmix(posterior, n_mcmc)

        cap        <- sum(hd[, 1] + hd[, 2])
        ccg_mean   <- ev_cur / n0
        postmean   <- unname(summary(posterior)["mean"])
        shrink_abs <- abs(ccg_mean - postmean)
        shrink_rel <- if(ccg_mean == 0) NA_real_ else shrink_abs / abs(ccg_mean)
        w_rob_post <- unname(posterior["w", ncol(posterior)])

        list(prior      = map,
             prior_rob  = map_robust,
             posterior  = posterior,
             pred_dist  = pred_dist,
             pred_int   = pred_int,
             eff_n_map  = .ess_safe(map,        cap),
             eff_n      = .ess_safe(map_robust, cap),
             eff_n_post = .ess_safe(posterior,  cap),
             shrinkage  = c(abs = shrink_abs, rel = shrink_rel),
             weights    = data.frame(w_prior = c(1 - robust, robust),
                                     w_post  = c(1 - w_rob_post, w_rob_post)),
             tau_prior  = .tau_prior_out(tau_dist, pr$tau_p),
             beta_prior = pr$beta_p,
             samples    = samples)
}

#-------------------------------------------------------------------------------
# Internal: SCIs of the treatment-vs-control contrasts from the posterior-sample
# matrix (binomial endpoint). Columns of 'simpost' (proportions) are in newdat
# row order. Returns list(dscir, cmatr).
#-------------------------------------------------------------------------------
.sci_binom <- function(simpost, nni, contr, type, base, alpha, alt_scs){
        if(contr == "RR"){
                cmatr <- mratios::contrMatRatio(n = nni, type = type, base = base)
                cc    <- BSagri::CCRatio.default(x = simpost, cmat = cmatr)
                SCIR  <- BSagri::SCSnp(x = cc, conf.level = 1 - alpha, alternative = alt_scs)
                est   <- SCIR$estimate; lo <- SCIR$conf.int[, 1]; up <- SCIR$conf.int[, 2]
        } else if(contr == "OR"){
                cmatr <- multcomp::contrMat(n = nni, type = type, base = base)
                cc    <- BSagri::CCDiff.default(x = log(simpost / (1 - simpost)), cmat = cmatr)
                SCIR  <- BSagri::SCSnp(x = cc, conf.level = 1 - alpha, alternative = alt_scs)
                est   <- exp(SCIR$estimate); lo <- exp(SCIR$conf.int[, 1]); up <- exp(SCIR$conf.int[, 2])
        } else {
                cmatr <- multcomp::contrMat(n = nni, type = type, base = base)
                cc    <- BSagri::CCDiff.default(x = simpost, cmat = cmatr)
                SCIR  <- BSagri::SCSnp(x = cc, conf.level = 1 - alpha, alternative = alt_scs)
                est   <- SCIR$estimate; lo <- SCIR$conf.int[, 1]; up <- SCIR$conf.int[, 2]
        }
        dscir <- data.frame(row.names = seq_along(SCIR$estimate),
                            comp     = names(SCIR$estimate),
                            estimate = unname(signif(est, 3)),
                            lower.ci = unname(signif(lo, 3)),
                            upper.ci = unname(signif(up, 3)))
        list(dscir = dscir, cmatr = cmatr)
}

#-------------------------------------------------------------------------------
# Internal: multi-arm SCI_binom().
#-------------------------------------------------------------------------------
.sci_binom_multiarm <- function(histdat, newdat, prior_est, contr, type, base,
                                robust, tau_prior, tau_dist, beta_prior,
                                limits, limits_pi, alpha, n_mcmc, seed){

        colnames(newdat)[1:3] <- c("group", "events", "non_events")
        arms        <- as.character(newdat$group)
        borrow_arms <- names(histdat)
        alt_scs     <- switch(limits, "lower" = "greater", "upper" = "less", "two.sided")

        cols   <- vector("list", length(arms))
        borrow <- vector("list", 0)
        for(i in seq_along(arms)){
                a <- arms[i]
                if(a %in% borrow_arms){
                        b <- .borrow_binom_arm(histdat[[a]], newdat$events[i], newdat$non_events[i],
                                               prior_est, robust, tau_prior, tau_dist, beta_prior,
                                               limits_pi, alpha, n_mcmc)
                        cols[[i]]   <- b$samples
                        borrow[[a]] <- b
                } else {
                        # Non-borrowing arm: Beta(1 + events, 1 + non-events).
                        cols[[i]] <- rbeta(n_mcmc, 1 + newdat$events[i], 1 + newdat$non_events[i])
                }
        }
        simpost <- do.call(cbind, cols); colnames(simpost) <- arms

        nni <- newdat$events + newdat$non_events; names(nni) <- newdat$group
        sci <- .sci_binom(simpost, nni, contr, type, base, alpha, alt_scs)

        newdat$mean <- newdat$events / (newdat$events + newdat$non_events)
        borrow <- borrow[intersect(arms, names(borrow))]
        pick   <- function(nm) lapply(borrow, `[[`, nm)
        pickv  <- function(nm) vapply(borrow, `[[`, numeric(1), nm)

        out_list <- list("comp"        = sci$dscir,
                         "contr"       = contr,
                         "histdat"     = histdat,
                         "newdat"      = newdat,
                         "distr"       = "betabinomial",
                         "borrow_arms" = names(borrow),
                         "prior_par"   = NA,
                         "weights"     = pick("weights"),
                         "eff_n_map"   = pickv("eff_n_map"),
                         "eff_n"       = pickv("eff_n"),
                         "eff_n_post"  = pickv("eff_n_post"),
                         "estimates"   = NA,
                         "contr_mat"   = sci$cmatr,
                         "prior"       = pick("prior"),
                         "prior_rob"   = pick("prior_rob"),
                         "posterior"   = pick("posterior"),
                         "pred_dist"   = pick("pred_dist"),
                         "pred_int"    = pick("pred_int"),
                         "shrinkage"   = pick("shrinkage"),
                         "simpost"     = simpost,
                         "tau_prior"   = pick("tau_prior"),
                         "beta_prior"  = pick("beta_prior"),
                         "alpha"       = alpha,
                         "seed"        = if(is.null(seed)) NA_integer_ else seed)

        structure(out_list, class = c("SCI", "MultiMAP_binom", "MultiMAP"))
}

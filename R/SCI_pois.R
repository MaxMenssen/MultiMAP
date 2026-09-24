#' Bayesian dynamic borrowing with simultaneous credible intervals (Poisson endpoint)
#'
#' Bayesian dynamic borrowing from historical control data for a Poisson (count)
#' endpoint with an exposure/offset. The pipeline: (1) build informative priors for the study arms of interest
#' based on external data (such as historical controls) -- either meta-analytic-predictive (MAP)
#' or an empirical-Bayes gamma prior; (2) robustify priors by mixing in a
#' non-informative component; (3) update to posteriors with the current
#' control data; (4) draw samples from the joint posterior; 
#' (5) derive SIMULTANEOUS credible sets for the comparisons of interest 
#' based on the joint posterior 
#'
#' @param histdat external (historical) data for the study-arms of interest; 
#' borrowing for one study arm (e.g. control group): submit a data.frame with columns
#'   \code{events} (the counted endpoint) and \code{offset} (the exposure / number of
#'   experimental units). 
#'   To borrow for more than one arm, pass a \emph{named list}
#'   containing different data frames: each element supplies the external
#'    data for a particular study-arm. 
#' @param newdat a data.frame containing the current trial; columns \code{group}, \code{events},
#'   \code{offset}. In case for borrowing for the control group: Row 1 is the concurrent 
#'   control, the remaining rows are treatment groups.
#' @param prior_est \code{"MAP"} (a meta-analytic-predictive prior via \code{RBesT::gMAP})
#'  or \code{"EB"} (an empirical-Bayes gamma prior following Tarone 1982)
#' @param type type of comparison: either \code{"Dunnett"} (many-to-one, each
#'   treatment vs. control) or \code{"Tukey"} (all pairwise comparisons).
#' @param contr \code{"rate_ratio"} provides ratios of rates (the default); 
#' \code{"rate_diff"} provides differences of rates.
#' @param limits limits of the credible intervals: either \code{"two.sided"},
#'   \code{"lower"} or \code{"upper"}.
#' @param limits_pi limits of the prediction interval for the concurrent control
#'   rate: either \code{"two.sided"}, \code{"lower"} or \code{"upper"}. Defaults
#'   to \code{limits} so that, unless set, the prediction interval follows the
#'   contrasts' one- or two-sidedness.
#' @param base index of the control group for the contrast matrix.
#' @param robust weight of the weakly-informative (robust) prior component.
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
#' @param alpha 1 - simultaneous confidence level.
#' @param n_mcmc total number of Monte-Carlo samples to draw.
#' @param seed optional integer; if supplied, \code{set.seed(seed)} is called so
#'   the \code{gMAP} MCMC (for \code{prior_est = "MAP"}) and the Monte-Carlo
#'   sampling are reproducible. \code{NULL} (default) leaves the RNG state
#'   untouched.
#' @param all_zero length-2 numeric giving the fallback gamma shape and rate
#'   \code{c(a, b)} used for the empirical-Bayes prior when all historical counts
#'   are zero (so the moment estimator is undefined).
#'
#' @return A MultiMAP object: a named list with S3 class
#'   \code{c("SCI", "MultiMAP_pois", "MultiMAP")} and the following slots:
#'   \describe{
#'     \item{\code{comp}}{data frame of the simultaneous credible intervals: one
#'       row per contrast with columns \code{comp} (label), \code{estimate},
#'       \code{lower.ci}, \code{upper.ci}.}
#'     \item{\code{contr}}{the contrast type used (\code{"rate_ratio"} or
#'       \code{"rate_diff"}).}
#'     \item{\code{histdat}}{the historical control data used.}
#'     \item{\code{newdat}}{the current trial data (concurrent control in row 1),
#'       with the observed per-group rate and treatment posterior gamma
#'       parameters added.}
#'     \item{\code{distr}}{the endpoint family, \code{"poisson"}.}
#'     \item{\code{prior_par}}{the gamma parameters (\code{a}, \code{b}) of the
#'       robustified prior mixture's components.}
#'     \item{\code{weights}}{data frame of the prior and posterior weights of the
#'       robustified mixture's components.}
#'     \item{\code{eff_n_map}}{effective sample size of the (non-robustified)
#'       prior.}
#'     \item{\code{eff_n}}{effective sample size of the robustified prior.}
#'     \item{\code{eff_n_post}}{effective sample size of the concurrent-control
#'       posterior.}
#'     \item{\code{estimates}}{the internal table of prior and posterior mixture
#'       parameters.}
#'     \item{\code{contr_mat}}{the contrast matrix used.}
#'     \item{\code{prior}}{the (non-robustified) prior (RBesT gamma mixture).}
#'     \item{\code{prior_rob}}{the robustified prior (RBesT gamma mixture).}
#'     \item{\code{posterior}}{the concurrent-control posterior (RBesT gamma
#'       mixture).}
#'     \item{\code{pred_dist}}{the prior-predictive distribution for the number of
#'       events in a future concurrent control group.}
#'     \item{\code{pred_int}}{the prediction interval for the concurrent control
#'       rate (named \code{median}, \code{lower}, \code{upper}).}
#'     \item{\code{shrinkage}}{named numeric \code{c(abs, rel)}: the absolute
#'       shrinkage of the concurrent control rate towards the borrowed prior and
#'       its relative version (absolute shrinkage divided by the raw concurrent
#'       control rate).}
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
#' @seealso \code{\link{SCI_norm}}, \code{\link{SCI_binom}}
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
#' Anscombe, F. J. (1950). Sampling theory of the negative binomial and
#' logarithmic series distributions. \emph{Biometrika, 37}(3/4), 358--382.
#' \doi{10.2307/2332388}
#'
#' Besag, J., Green, P., Higdon, D., & Mengersen, K. (1995). Bayesian computation
#' and stochastic systems. \emph{Statistical Science, 10}(1), 3--66.
#' \url{http://www.jstor.org/stable/2246224}
#'
#' @examples
#' \donttest{
#' # Egg-laying count data from EFSA (2025), Annex C: number of eggs laid by a
#' # control group of 12 hens over 10 weeks. HCD = 106 historical control studies
#' # summarised by their mean eggs/hen (Table 3); index study = total eggs per
#' # dose group (Table 2). Poisson mapping: offset = 12 hens,
#' # events = total eggs = round(mean eggs/hen * 12), so rate = eggs/hen.
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
#' # Index study (Table 2): total eggs per dose group; dose 0 = concurrent control
#' newdat <- data.frame(group  = c("0", "100", "200", "400"),  # mg/day
#'                      events = c(546, 419, 473, 420),
#'                      offset = 12)
#'
#' # Simultaneous credible intervals for the dose-vs-control egg-laying rate ratios
#' res <- SCI_pois(histdat, newdat, prior_est = "MAP",
#'                 contr = "rate_ratio", seed = 81771)
#' summary(res)
#' }
#' @export
SCI_pois <- function(histdat,
                     newdat,
                     prior_est = "MAP",
                     type = "Dunnett",
                     contr = "rate_ratio",
                     limits = "two.sided",
                     limits_pi = limits,
                     base = 1,
                     robust = 0.2,
                     tau_prior = NULL,
                     tau_dist = "HalfNormal",
                     beta_prior = NULL,
                     alpha = 0.05,
                     n_mcmc = 50000,
                     seed = NULL,
                     all_zero = c(0.0001, 1)){

        # Optional reproducibility: seeds gMAP()'s MCMC (MAP path) and the
        # Monte-Carlo sampling below.
        if(!is.null(seed)) set.seed(seed)

        # User-facing one-sidedness vocabulary: two.sided / lower / upper.
        # 'limits' drives the contrasts, 'limits_pi' the prediction
        # interval (independent; defaults to 'limits').
        limits    <- match.arg(limits,    c("two.sided", "lower", "upper"))
        limits_pi <- match.arg(limits_pi, c("two.sided", "lower", "upper"))

        # Contrast type: ratio of rates / difference of rates.
        contr <- match.arg(contr, c("rate_ratio", "rate_diff"))

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
                stop("'histdat' must be a data frame with at least 2 columns (events, offset).")
        if(!is.data.frame(newdat) || ncol(newdat) < 3)
                stop("'newdat' must be a data frame with at least 3 columns (group, events, offset).")
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
        # remaining arms use the non-borrowing near-flat Gamma(eps, eps) posterior.
        if(multi_arm){
                .check_multi_histdat(histdat, newdat[[1]])
                return(.sci_pois_multiarm(histdat, newdat, prior_est, contr, type, base,
                                          robust, tau_prior, tau_dist, beta_prior,
                                          limits, limits_pi, alpha, n_mcmc, all_zero, seed))
        }

        # Near-flat gamma prior for the (non-borrowing) treatment arms.
        eps <- 0.001

        # Effective gMAP hyperpriors actually used (filled in the MAP block; NA
        # for the EB path). Stored in the output.
        tau_p  <- NA
        beta_p <- NA

        #=======================================================================
        # Prior specification: MAP (gMAP) or EB (empirical Bayes, method of moments)
        #=======================================================================

        if(prior_est == "MAP"){

                histdat_map <- cbind(data.frame(study = factor(seq_len(nrow(histdat)))),
                                     histdat)
                colnames(histdat_map)[2:3] <- c("events", "offset")

                # Pooled historical (grand) rate; its log centres the intercept prior.
                lambda_bar <- mean(histdat_map$events / histdat_map$offset)

                if(!is.finite(lambda_bar) || lambda_bar <= 0)
                        stop("The MAP prior needs positive historical counts; use prior_est = 'EB' for all-zero data.")

                # Default gMAP hyperpriors on the log-rate scale (both overridable):
                #  - tau ~ HalfNormal(0, 0.5): a fixed, weakly-informative between-study
                #    heterogeneity prior. On the log-rate scale the heterogeneity is a
                #    coefficient of variation and is (to first order) invariant to the
                #    rate, so a fixed scale is preferable to one tied to the rate/count
                #    (Roever et al. 2021). cbind(0, 1) gives a more diffuse alternative.
                #  - intercept ~ Normal(log(lambda_bar), 1): centred at the pooled
                #    historical log-rate with a weakly-informative SD (use 2 for a more
                #    diffuse prior). Centring at 0 (rate 1) would badly miscalibrate the
                #    MAP prior for rates far from 1.
                tau_p  <- if(is.null(tau_prior))  cbind(0, 0.5)             else tau_prior
                beta_p <- if(is.null(beta_prior)) cbind(log(lambda_bar), 1) else beta_prior

                map_mcmc <- gMAP(events ~ 1 + offset(log(offset)) | study,
                                 data       = histdat_map,
                                 family     = poisson,
                                 tau.dist   = tau_dist,
                                 tau.prior  = tau_p,
                                 beta.prior = beta_p)

                # Derive gamma mixture on the rate scale
                map <- automixfit(map_mcmc, type = "gamma")
        }

        if(prior_est == "EB"){

                # No historical events -> fall back to the supplied gamma prior.
                if(all(histdat[, 1] == 0)){
                        a_est <- all_zero[1]
                        b_est <- all_zero[2]
                }
                else{
                        # Empirical-Bayes gamma prior following Tarone (1982): the
                        # control rate varies across historical studies as
                        # Gamma(shape = p, rate = p / mu), so borrowing the
                        # historical data into the concurrent control is the
                        # conjugate gamma-Poisson update that adds p prior events
                        # over p / mu prior exposure units. The prior mean mu is the
                        # pooled historical rate; the negative-binomial shape p
                        # (which sets the borrowing strength) is estimated here by
                        # the method of moments -- a closed-form, convergence-free
                        # alternative to the negative-binomial ML fit Tarone used
                        # (glm.nb is unreliable on small and/or underdispersed data).
                        #
                        # Homogeneous historical counts -> large p -> strong
                        # borrowing; heterogeneous counts -> small p -> little
                        # borrowing. Following Tarone's recommendation, when the
                        # data are equi-/under-dispersed (the moment denominator is
                        # <= 0) or the implied prior exposure p / mu would exceed the
                        # total historical exposure, the prior is capped at the
                        # pooled historical data Gamma(sum(events), sum(offset)),
                        # i.e. full pooling -- the most informative admissible prior.
                        # (This replaces the diffuse prior a naive moment floor
                        # returns for under-dispersion, and mirrors the constrained
                        # moment estimator used for the binomial endpoint.)
                        y_h    <- histdat[, 1]
                        off_h  <- histdat[, 2]
                        lambda_hat <- sum(y_h) / sum(off_h)        # pooled rate (mu)
                        mu_h       <- lambda_hat * off_h           # expected counts
                        denom      <- sum((y_h - mu_h)^2 - mu_h)   # excess-variance sum

                        if(denom <= 0){
                                # equi-/under-dispersed -> full pooling (Tarone cap)
                                a_est <- sum(y_h)
                                b_est <- sum(off_h)
                                warning("HCD is equi-/under-dispersed: the prior was capped at the pooled historical data (full pooling).")
                        }
                        else{
                                p_hat <- sum(mu_h^2) / denom       # negative-binomial shape
                                a_est <- p_hat                     # prior events (shape)
                                b_est <- p_hat / lambda_hat        # prior exposure (rate = p / mu)

                                # Tarone: p / mu can exceed the total exposure for very
                                # homogeneous data -> cap at the pooled historical data.
                                if(b_est > sum(off_h)){
                                        a_est <- sum(y_h)
                                        b_est <- sum(off_h)
                                }
                        }
                }

                map <- mixgamma(c(1, a_est, b_est), param = "ab", likelihood = "poisson")
        }

        #=======================================================================
        # Robustify, update, sample
        #=======================================================================

        # Prior mean rate (used to centre the weakly-informative component)
        prior_mean <- unname(summary(map)["mean"])

        map_robust <- robustify(map, weight = robust, mean = prior_mean)

        prior_tab <- data.frame(cbind(t(map_robust)))
        colnames(prior_tab) <- c("w", "a", "b")

        #-----------------------------------------------------------------------
        # Prediction interval for the concurrent control rate.
        # preddist() on the (non-robustified) prior gives the gamma-Poisson
        # predictive for the number of events over the control exposure; dividing
        # the count quantiles by that exposure puts the interval on the rate scale.

        offset0   <- newdat[1, 3]
        pred_dist <- RBesT::preddist(map, n = offset0)

        if(limits_pi == "two.sided"){
                pred_int <- c(median = unname(RBesT::qmix(pred_dist, 0.5))        / offset0,
                              lower  = unname(RBesT::qmix(pred_dist, alpha/2))     / offset0,
                              upper  = unname(RBesT::qmix(pred_dist, 1 - alpha/2)) / offset0)
        }

        if(limits_pi == "lower"){
                pred_int <- c(median = unname(RBesT::qmix(pred_dist, 0.5))   / offset0,
                              lower  = unname(RBesT::qmix(pred_dist, alpha))  / offset0,
                              upper  = NA)
        }

        if(limits_pi == "upper"){
                pred_int <- c(median = unname(RBesT::qmix(pred_dist, 0.5))        / offset0,
                              lower  = NA,
                              upper  = unname(RBesT::qmix(pred_dist, 1 - alpha))  / offset0)
        }

        # Translate 'limits' to the BSagri::SCSnp() vocabulary.
        alt_scs <- switch(limits,
                          "lower" = "greater",
                          "upper" = "less",
                          "two.sided")

        #-----------------------------------------------------------------------
        # Posterior for the concurrent control rate: update the robustified prior
        # with the control count over its exposure (postmix's n = exposure,
        # m = count / exposure, so a + count and b + exposure are added).

        posterior <- postmix(map_robust,
                             n = offset0,
                             m = newdat[1, 2] / offset0)

        post_tab <- data.frame(cbind(t(posterior)))
        colnames(post_tab) <- c("w_post", "post_a", "post_b")

        #-----------------------------------------------------------------------
        # Sample the control posterior mixture (rates)

        simpost_ccg <- RBesT::rmix(posterior, n_mcmc)

        #-----------------------------------------------------------------------
        # Treatment posteriors: near-flat Gamma(eps, eps) prior updated with the
        # arm's count and exposure -> Gamma(eps + events, eps + offset); no
        # historical information is borrowed for the treatment arms.

        newdat$posta <- eps + newdat[, 2]
        newdat$postb <- eps + newdat[, 3]

        simpost_trt <- apply(X      = newdat[-1, c("posta", "postb")],
                             MARGIN = 1,
                             FUN    = function(x){rgamma(n     = length(simpost_ccg),
                                                         shape = x["posta"],
                                                         rate  = x["postb"])})

        #-----------------------------------------------------------------------
        # Matrix with samples from posteriors (control first)

        simpost <- cbind(simpost_ccg,
                         simpost_trt)

        #-----------------------------------------------------------------------
        # Output object for weights (robust component is the last mixture component)

        w_rob_post <- post_tab$w_post[nrow(post_tab)]
        weights <- data.frame(w_prior = c(1 - robust, robust),
                              w_post  = c(1 - w_rob_post, w_rob_post))

        #-----------------------------------------------------------------------
        # Prior parameters + estimates output

        a_est <- prior_tab$a
        b_est <- prior_tab$b

        estimates_out <- cbind(prior_tab, post_tab)

        # Observed per-group rate (concurrent control in row 1)
        newdat$mean <- newdat[, 2] / newdat[, 3]

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
        # Effective sample size of the concurrent-control posterior

        eff_n_post <- try(round(ess(posterior, method = "elir")), silent = TRUE)

        if(inherits(eff_n_post, "try-error")){

                eff_n_post <- round(ess(posterior, method = "moment"))

                warning("Estimation of ESS via elir-method does not work.\n ESS was estimated based on method of moments")
        }

        eff_n_post <- min(max(0, eff_n_post), sum(histdat[, 2]))

        #-----------------------------------------------------------------------
        # Shrinkage

        postmean <- unname(summary(posterior)["mean"])
        ccg_mean <- newdat[1, 2] / newdat[1, 3]

        shrink_abs <- unname(abs(ccg_mean - postmean))
        shrink_rel <- if(ccg_mean == 0) NA_real_ else shrink_abs / abs(ccg_mean)
        shrinkage  <- c(abs = shrink_abs, rel = shrink_rel)

        #-----------------------------------------------------------------------
        ## ratio of rates (rate_1 / rate_2)

        if(contr == "rate_ratio"){
                nni <- newdat[, 3]
                names(nni) <- newdat[, 1]

                cmatr <- mratios::contrMatRatio(n    = nni,
                                                type = type,
                                                base = base)

                ratio2control <- BSagri::CCRatio.default(x = simpost, cmat = cmatr)

                SCIR <- BSagri::SCSnp(x = ratio2control, conf.level = 1 - alpha,
                                      alternative = alt_scs)

                dscir <- data.frame(
                        row.names = seq_along(SCIR$estimate),
                        comp      = names(SCIR$estimate),
                        estimate  = unname(signif(SCIR$estimate, 3)),
                        lower.ci  = unname(signif(SCIR$conf.int[, 1], 3)),
                        upper.ci  = unname(signif(SCIR$conf.int[, 2], 3)))
        }

        #-----------------------------------------------------------------------
        ## difference of rates (rate_1 - rate_2)

        if(contr == "rate_diff"){
                nni <- newdat[, 3]
                names(nni) <- newdat[, 1]

                cmatr <- multcomp::contrMat(n    = nni,
                                            type = type,
                                            base = base)

                diff2control <- BSagri::CCDiff.default(x = simpost, cmat = cmatr)

                SCIR <- BSagri::SCSnp(x = diff2control, conf.level = 1 - alpha,
                                      alternative = alt_scs)

                dscir <- data.frame(
                        row.names = seq_along(SCIR$estimate),
                        comp      = names(SCIR$estimate),
                        estimate  = unname(signif(SCIR$estimate, 3)),
                        lower.ci  = unname(signif(SCIR$conf.int[, 1], 3)),
                        upper.ci  = unname(signif(SCIR$conf.int[, 2], 3)))
        }

        #-----------------------------------------------------------------------
        # Output

        out_list <- list("comp"       = dscir,
                         "contr"      = contr,
                         "histdat"    = histdat,
                         "newdat"     = newdat,
                         "distr"      = "poisson",
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
                            class = c("SCI", "MultiMAP_pois", "MultiMAP"))

        return(out_s3)
}


#-------------------------------------------------------------------------------
# Internal: build the informative prior (RBesT gamma mixture) for one arm from its
# historical data (events, offset), via MAP (gMAP) or EB (empirical Bayes, Tarone
# 1982). Shared by SCI_pois() and PI_pois(). Returns list(map, tau_p, beta_p).
#-------------------------------------------------------------------------------
.map_prior_pois <- function(hd, prior_est, tau_prior, tau_dist, beta_prior, all_zero){

        tau_p <- NA; beta_p <- NA

        if(prior_est == "MAP"){
                hd_map <- cbind(data.frame(study = factor(seq_len(nrow(hd)))), hd)
                colnames(hd_map)[2:3] <- c("events", "offset")

                lambda_bar <- mean(hd_map$events / hd_map$offset)
                if(!is.finite(lambda_bar) || lambda_bar <= 0)
                        stop("The MAP prior needs positive historical counts; use prior_est = 'EB' for all-zero data.")

                tau_p  <- if(is.null(tau_prior))  cbind(0, 0.5)             else tau_prior
                beta_p <- if(is.null(beta_prior)) cbind(log(lambda_bar), 1) else beta_prior

                map_mcmc <- gMAP(events ~ 1 + offset(log(offset)) | study,
                                 data       = hd_map,
                                 family     = poisson,
                                 tau.dist   = tau_dist,
                                 tau.prior  = tau_p,
                                 beta.prior = beta_p)
                map <- automixfit(map_mcmc, type = "gamma")
        } else {
                # Empirical-Bayes gamma prior following Tarone (1982); see SCI_pois().
                if(all(hd[, 1] == 0)){
                        a_est <- all_zero[1]
                        b_est <- all_zero[2]
                } else {
                        y_h    <- hd[, 1]
                        off_h  <- hd[, 2]
                        lambda_hat <- sum(y_h) / sum(off_h)
                        mu_h       <- lambda_hat * off_h
                        denom      <- sum((y_h - mu_h)^2 - mu_h)

                        if(denom <= 0){
                                a_est <- sum(y_h)
                                b_est <- sum(off_h)
                                warning("HCD is equi-/under-dispersed: the prior was capped at the pooled historical data (full pooling).")
                        } else {
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
        list(map = map, tau_p = tau_p, beta_p = beta_p)
}

#-------------------------------------------------------------------------------
# Internal: borrowing pipeline for ONE arm (Poisson endpoint). (ev_cur, off_cur)
# is the arm's current-trial count and exposure. Returns the per-arm bundle with
# n_mcmc posterior draws (rates) in $samples.
#-------------------------------------------------------------------------------
.borrow_pois_arm <- function(hd, ev_cur, off_cur, prior_est, robust,
                             tau_prior, tau_dist, beta_prior, limits_pi, alpha, n_mcmc, all_zero){

        pr  <- .map_prior_pois(hd, prior_est, tau_prior, tau_dist, beta_prior, all_zero)
        map <- pr$map
        prior_mean <- unname(summary(map)["mean"])
        map_robust <- robustify(map, weight = robust, mean = prior_mean)

        pred_dist <- RBesT::preddist(map, n = off_cur)
        pred_int  <- .pred_int_from(pred_dist, limits_pi, alpha, scale = off_cur)

        posterior <- postmix(map_robust, n = off_cur, m = ev_cur / off_cur)
        samples   <- RBesT::rmix(posterior, n_mcmc)

        cap        <- sum(hd[, 2])
        ccg_mean   <- ev_cur / off_cur
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
# matrix (Poisson endpoint). Columns of 'simpost' (rates) are in newdat row order.
# Returns list(dscir, cmatr).
#-------------------------------------------------------------------------------
.sci_pois <- function(simpost, nni, contr, type, base, alpha, alt_scs){
        if(contr == "rate_ratio"){
                cmatr <- mratios::contrMatRatio(n = nni, type = type, base = base)
                cc    <- BSagri::CCRatio.default(x = simpost, cmat = cmatr)
        } else {
                cmatr <- multcomp::contrMat(n = nni, type = type, base = base)
                cc    <- BSagri::CCDiff.default(x = simpost, cmat = cmatr)
        }
        SCIR <- BSagri::SCSnp(x = cc, conf.level = 1 - alpha, alternative = alt_scs)
        dscir <- data.frame(row.names = seq_along(SCIR$estimate),
                            comp     = names(SCIR$estimate),
                            estimate = unname(signif(SCIR$estimate, 3)),
                            lower.ci = unname(signif(SCIR$conf.int[, 1], 3)),
                            upper.ci = unname(signif(SCIR$conf.int[, 2], 3)))
        list(dscir = dscir, cmatr = cmatr)
}

#-------------------------------------------------------------------------------
# Internal: multi-arm SCI_pois().
#-------------------------------------------------------------------------------
.sci_pois_multiarm <- function(histdat, newdat, prior_est, contr, type, base,
                               robust, tau_prior, tau_dist, beta_prior,
                               limits, limits_pi, alpha, n_mcmc, all_zero, seed){

        colnames(newdat)[1:3] <- c("group", "events", "offset")
        arms        <- as.character(newdat$group)
        borrow_arms <- names(histdat)
        alt_scs     <- switch(limits, "lower" = "greater", "upper" = "less", "two.sided")

        # Near-flat gamma prior for the non-borrowing arms.
        eps <- 0.001

        cols   <- vector("list", length(arms))
        borrow <- vector("list", 0)
        for(i in seq_along(arms)){
                a <- arms[i]
                if(a %in% borrow_arms){
                        b <- .borrow_pois_arm(histdat[[a]], newdat$events[i], newdat$offset[i],
                                              prior_est, robust, tau_prior, tau_dist, beta_prior,
                                              limits_pi, alpha, n_mcmc, all_zero)
                        cols[[i]]   <- b$samples
                        borrow[[a]] <- b
                } else {
                        # Non-borrowing arm: Gamma(eps + events, eps + offset).
                        cols[[i]] <- rgamma(n_mcmc, shape = eps + newdat$events[i],
                                            rate = eps + newdat$offset[i])
                }
        }
        simpost <- do.call(cbind, cols); colnames(simpost) <- arms

        nni <- newdat$offset; names(nni) <- newdat$group
        sci <- .sci_pois(simpost, nni, contr, type, base, alpha, alt_scs)

        newdat$mean <- newdat$events / newdat$offset
        borrow <- borrow[intersect(arms, names(borrow))]
        pick   <- function(nm) lapply(borrow, `[[`, nm)
        pickv  <- function(nm) vapply(borrow, `[[`, numeric(1), nm)

        out_list <- list("comp"        = sci$dscir,
                         "contr"       = contr,
                         "histdat"     = histdat,
                         "newdat"      = newdat,
                         "distr"       = "poisson",
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

        structure(out_list, class = c("SCI", "MultiMAP_pois", "MultiMAP"))
}

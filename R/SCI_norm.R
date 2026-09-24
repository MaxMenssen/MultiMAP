#' Bayesian dynamic borrowing with simultaneous credible intervals (normal endpoint)
#'
#' Bayesian dynamic borrowing from historical control data for a normally
#' distributed endpoint. The pipeline: (1) build informative meta-analytic-predictive 
#' (MAP) priors for the study arms of interest based on external data (such as 
#' historical controls); (2) robustify priors by mixing in a
#' non-informative component; (3) update to posteriors with the current
#' control data; (4) draw samples from the joint posterior; 
#' (5) derive SIMULTANEOUS credible sets for the comparisons of interest 
#' based on the joint posterior 
#' 
#'
#' @param histdat external (historical) data for the study-arms of interest; 
#' borrowing for one study arm (e.g. control group): data.frame with one row per study, 
#' columns \code{mean}, \code{sd}, \code{n}. To borrow for more than one arm, 
#' pass a \emph{named list} containing different data frames: each element supplies 
#' the external data for a particular study-arm. 
#' @param newdat a data.frame containing the current trial; columns \code{group}, \code{mean}, \code{sd},
#'   \code{n}. In case for borrowing for the control group: Row 1 is the concurrent 
#'   control, the remaining rows are treatment groups.
#' @param prior_est currently, only \code{"MAP"} is supported: a
#'   meta-analytic-predictive prior via \code{RBesT::gMAP} (\code{family =
#'   gaussian}), fitted to a normal mixture.
#' @param type type of comparison: either \code{"Dunnett"} (many-to-one, each
#'   treatment vs. control) or \code{"Tukey"} (all pairwise comparisons).
#' @param contr \code{"mean_diff"} gives treatment-minus-control mean
#'   differences; \code{"mean_ratio"} gives ratios of means (requires positive
#'   means).
#' @param limits limits of the credible intervals: either \code{"two.sided"},
#'   \code{"lower"} or \code{"upper"}.
#' @param limits_pi limits of the prediction interval for the concurrent control
#'   mean: either \code{"two.sided"}, \code{"lower"} or \code{"upper"}. Default is
#'   set to \code{limits} so, unless set, the prediction interval follows the
#'   contrasts' one- or two-sidedness.
#' @param base index of the control group for the contrast matrix.
#' @param robust weight of the weakly-informative (robust) prior component.
#' @param sigma reference scale for the RBesT Gaussian machinery. If \code{NULL},
#'   \code{sigma} is estimated from the data as an n-weighted pooled sd of the
#'   historical controls (unit-information standard deviation, Röver et al. 2021); a supplied value 
#'   is used as a fixed known reference scale throughout.
#' @param tau_prior optional override for the \code{gMAP} between-study
#'   heterogeneity prior (\code{tau.prior} in \code{RBesT::gMAP}). \code{NULL}
#'   (the default) uses the built-in heuristic (a \code{HalfNormal(0, sigma/2)}).
#'   Supply any \code{cbind(0, scale)} to override.
#' @param tau_dist distribution family of the \code{gMAP} between-study
#'   heterogeneity prior (\code{tau.dist} in \code{RBesT::gMAP}); one of
#'   \code{"HalfNormal"} (the default), \code{"TruncNormal"}, \code{"Uniform"},
#'   \code{"Gamma"}, \code{"InvGamma"}, \code{"LogNormal"}, \code{"TruncCauchy"},
#'   \code{"Exp"} or \code{"Fixed"}. The built-in default \code{tau_prior}
#'   heuristic is calibrated for the half-normal, so when \code{tau_dist} is set
#'   to any other family a matching \code{tau_prior} must also be supplied.
#' @param beta_prior optional override for the \code{gMAP} intercept prior
#'   (\code{beta.prior} in \code{RBesT::gMAP}). \code{NULL} (the default) uses the
#'   built-in heuristic (\code{cbind(grand mean, sigma)}); supply e.g.
#'   \code{cbind(mean, sd)}.
#' @param alpha 1 - simultaneous confidence level.
#' @param n_mcmc total number of Monte-Carlo samples to draw.
#' @param seed optional integer; if supplied, \code{set.seed(seed)} is called so
#'   the \code{gMAP} MCMC and the Monte-Carlo sampling are reproducible.
#'   \code{NULL} (default) leaves the RNG state untouched.
#'
#' @return A MultiMAP object: a named list with S3 class
#'   \code{c("SCI", "MultiMAP_norm", "MultiMAP")} and the following slots:
#'   \describe{
#'     \item{\code{comp}}{data frame of the simultaneous credible intervals: one
#'       row per contrast with columns \code{comp} (label), \code{estimate},
#'       \code{lower.ci}, \code{upper.ci}.}
#'     \item{\code{contr}}{the contrast type used (\code{"mean_diff"} or
#'       \code{"mean_ratio"}).}
#'     \item{\code{histdat}}{the historical control data used (mean, sd, n).}
#'     \item{\code{newdat}}{the current trial data (concurrent control in row 1),
#'       with the per-group posterior mean/sd added.}
#'     \item{\code{distr}}{the endpoint family, \code{"normal"}.}
#'     \item{\code{prior_par}}{the sampled prior mixture parameters (means and
#'       sds of the components).}
#'     \item{\code{weights}}{data frame of the prior and posterior weights of the
#'       robustified mixture's components.}
#'     \item{\code{eff_n_map}}{effective sample size of the (non-robustified) MAP
#'       prior.}
#'     \item{\code{eff_n}}{effective sample size of the robustified MAP prior.}
#'     \item{\code{eff_n_post}}{effective sample size of the concurrent-control
#'       posterior.}
#'     \item{\code{estimates}}{the internal sampling table (mixture weights,
#'       means, sds and per-component draw counts).}
#'     \item{\code{contr_mat}}{the contrast matrix used.}
#'     \item{\code{prior}}{the MAP prior (RBesT normal mixture).}
#'     \item{\code{prior_rob}}{the robustified MAP prior (RBesT normal mixture).}
#'     \item{\code{posterior}}{the concurrent-control posterior (RBesT normal
#'       mixture).}
#'     \item{\code{pred_dist}}{the MAP-based predictive distribution for a future
#'       concurrent control mean.}
#'     \item{\code{pred_int}}{the prediction interval (named \code{median},
#'       \code{lower}, \code{upper}).}
#'     \item{\code{shrinkage}}{named numeric \code{c(abs, rel)}: the absolute
#'       shrinkage of the concurrent control mean towards the borrowed prior and
#'       its relative version (absolute shrinkage divided by the raw concurrent
#'       control mean).}
#'     \item{\code{simpost}}{matrix of posterior Monte-Carlo samples (control plus
#'       one column per treatment group).}
#'     \item{\code{tau_prior}}{a list with the gMAP between-study heterogeneity
#'       prior distribution name (\code{dist}) and its parameters (\code{prior})
#'       actually used (the defaults or the supplied overrides).}
#'     \item{\code{beta_prior}}{the gMAP intercept prior actually used (the
#'       default or the supplied override).}
#'     \item{\code{alpha}}{the \code{alpha} used.}
#'     \item{\code{seed}}{the \code{seed} used (\code{NA} if none).}
#'   }
#'   Use \code{\link{plot.MultiMAP}}, \code{\link{print.MultiMAP}},
#'   \code{\link{summary.MultiMAP}}, \code{\link{get_output}} and
#'   \code{\link{report_norm}} on it.
#'
#' @seealso \code{\link{PI_norm}}, \code{\link{SIM_norm}}, \code{\link{report_norm}}
#'
#' @references
#' EFSA (European Food Safety Authority) (2025). Use and reporting of historical
#' control data for regulatory studies. \emph{EFSA Journal}.
#' \doi{10.2903/j.efsa.2025.9576}
#'
#' Roever, C., Bender, R., Dias, S., Schmid, C. H., Schmidli, H., Sturtz, S.,
#' Weber, S., & Friede, T. (2021). On weakly informative prior distributions for
#' the heterogeneity parameter in Bayesian random-effects meta-analysis.
#' \emph{Research Synthesis Methods, 12}(4), 448--474. \doi{10.1002/jrsm.1475}
#' 
#' Schmidli, H., Gsteiger, S., Roychoudhury, S., O'Hagan, A., Spiegelhalter, D.,
#' & Neuenschwander, B. (2014). Robust meta-analytic-predictive priors in
#' clinical trials with historical control information. \emph{Biometrics,
#' 70}(4), 1023--1032. \doi{10.1111/biom.12242}
#'
#' Besag, J., Green, P., Higdon, D., & Mengersen, K. (1995). Bayesian computation
#' and stochastic systems. \emph{Statistical Science, 10}(1), 3--66.
#' \url{http://www.jstor.org/stable/2246224}
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
#' # SCI_norm() needs the per-study SD: sd = se * sqrt(n)
#' histdat <- data.frame(mean = hcd$mean, sd = hcd$se * sqrt(hcd$n), n = hcd$n)
#'
#' # Index study (EFSA 2025, Table 14): concurrent control (0 ppm) and three
#' # dose groups (ppm); row 1 is the concurrent control.
#' idx <- data.frame(group = c("0", "600", "1400", "5000"),
#'                   mean  = c(1.05, 1.07, 1.37, 1.84),
#'                   se    = c(0.138, 0.140, 0.164, 0.204),
#'                   n     = c(5, 5, 5, 5))
#' newdat <- data.frame(group = idx$group, mean = idx$mean,
#'                      sd = idx$se * sqrt(idx$n), n = idx$n)
#'
#' # Simultaneous credible intervals (treatment-vs-control mean differences)
#' res <- SCI_norm(histdat, newdat, contr = "mean_diff", seed = 81771)
#' summary(res)
#' }
#' @export
SCI_norm <- function(histdat,
                          newdat,
                          prior_est = "MAP",
                          type = "Dunnett",
                          contr = "mean_diff",
                          limits = "two.sided",
                          limits_pi = limits,
                          base = 1,
                          robust = 0.2,
                          sigma = NULL,
                          tau_prior = NULL,
                          tau_dist = "HalfNormal",
                          beta_prior = NULL,
                          alpha = 0.05,
                          n_mcmc = 50000,
                          seed = NULL){

        # Optional reproducibility: seeds gMAP()'s MCMC and the Monte-Carlo
        # sampling below.
        if(!is.null(seed)) set.seed(seed)

        # User-facing one-sidedness vocabulary: two.sided / lower / upper.
        # 'limits' drives the contrasts, 'limits_pi' the prediction
        # interval (independent; defaults to 'limits').
        limits    <- match.arg(limits,    c("two.sided", "lower", "upper"))
        limits_pi <- match.arg(limits_pi, c("two.sided", "lower", "upper"))

        # Contrast type: mean difference / mean ratio.
        contr <- match.arg(contr, c("mean_diff", "mean_ratio"))

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

        if(!is.data.frame(newdat) || ncol(newdat) < 4)
                stop("'newdat' must be a data frame with at least 4 columns (group, mean, sd, n).")
        if(!multi_arm && (!is.data.frame(histdat) || ncol(histdat) < 3))
                stop("'histdat' must be a data frame with at least 3 columns (mean, sd, n).")
        if(!is.numeric(robust) || length(robust) != 1 || robust < 0 || robust > 1)
                stop("'robust' must be a single number in [0, 1].")
        if(!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1)
                stop("'alpha' must be a single number in (0, 1).")
        if(!is.numeric(base) || length(base) != 1 || base < 1 || base > nrow(newdat))
                stop("'base' must be a single group index in 1:nrow(newdat).")

        #-----------------------------------------------------------------------
        # Column bookkeeping (mean, sd, n)

        colnames(newdat)[1:4]  <- c("group", "mean", "sd", "n")
        if(!multi_arm) colnames(histdat)[1:3] <- c("mean", "sd", "n")

        # Only the MAP prior is supported (the pooled approach has been removed).
        if(prior_est != "MAP"){
                stop("SCI_norm() supports only prior_est = 'MAP'.")
        }

        #-----------------------------------------------------------------------
        # Multi-arm borrowing: each named arm gets its own MAP prior + posterior;
        # the remaining arms use the non-borrowing (flat-prior) posterior. Handled
        # in a dedicated path so the single-arm code below is unchanged.
        if(multi_arm){
                .check_multi_histdat(histdat, newdat$group)
                return(.sci_norm_multiarm(histdat, newdat, contr, type, base,
                                          robust, sigma, tau_prior, tau_dist, beta_prior,
                                          limits, limits_pi, alpha, n_mcmc, seed))
        }

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

        # Effective gMAP hyperpriors actually used (filled in the MAP block below;
        # stored in the output for reproducibility).
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
                # Default tau ~ HalfNormal(0, sigma/2) follows the n_infinity
                # heuristic; the default beta.prior is the RBesT-recommended
                # unit-information prior on the intercept (SD = sigma), centred at
                # the historical grand mean. Both can be overridden via tau_prior /
                # beta_prior.

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
        
        pred_dist <- RBesT::preddist(map, n = newdat$n[1], sigma = sigma_ref)   # n = concurrent control sample size

        if(limits_pi == "two.sided"){
                pred_int <- c(median = unname(RBesT::qmix(pred_dist, 0.5)),
                      lower  = unname(RBesT::qmix(pred_dist, alpha/2)),
                      upper  = unname(RBesT::qmix(pred_dist, 1 - alpha/2)))

        }

        if(limits_pi == "lower"){
                pred_int <- c(median = unname(RBesT::qmix(pred_dist, 0.5)),
                      lower  = unname(RBesT::qmix(pred_dist, alpha)),
                      upper  = NA)

        }

        if(limits_pi == "upper"){
                pred_int <- c(median = unname(RBesT::qmix(pred_dist, 0.5)),
                      lower  = NA,
                      upper  = unname(RBesT::qmix(pred_dist, 1 - alpha)))

        }

        # Translate the user-facing 'limits' (two.sided / lower / upper) to
        # the BSagri::SCSnp() vocabulary (two.sided / greater / less): a lower
        # simultaneous bound corresponds to 'greater', an upper bound to 'less'.
        alt_scs <- switch(limits,
                          "lower" = "greater",
                          "upper" = "less",
                          "two.sided")
        
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
        # Posterior for the concurrent control group

        ccg_se   <- newdat$sd[1] / sqrt(newdat$n[1])
        posterior <- postmix(map_robust, m = newdat$mean[1], se = ccg_se)

        post_df <- data.frame(cbind(t(posterior)))
        colnames(post_df) <- c("w_post", "post_m", "post_s")
        post_df <- post_df[, c("post_m", "post_s", "w_post")]

        simdat <- data.frame(cbind(map_robust_df, post_df))
        colnames(simdat)[1:3] <- c("w", "m", "s")
        simdat$n       <- n_mcmc                          # total no. of samples to draw
        simdat$n_prior <- round(simdat$w      * simdat$n) # samples from each prior comp.
        simdat$n_post  <- round(simdat$w_post * simdat$n) # samples from each posterior comp.


        #-----------------------------------------------------------------------
        # Sample from the normal posterior mixture (control group)

        simpost_ccg_fun <- function(x){
                rnorm(n    = x$n_post,
                      mean = x$post_m,
                      sd   = x$post_s)
        }

        simpost_ccg <- unlist(lapply(X   = split(simdat, seq(nrow(simdat))),
                                     FUN = simpost_ccg_fun))


        #-----------------------------------------------------------------------
        # Posterior for the treatment groups
        # with improper prior on the mean (all means are equally likely)):
        # mean_g ~ N(ybar_g, sd_g / sqrt(n_g))

        newdat$post_m <- newdat$mean
        newdat$post_s <- newdat$sd / sqrt(newdat$n)

        simpost_trt <- apply(X      = newdat[-1, c("post_m", "post_s")],
                             MARGIN = 1,
                             FUN    = function(x){rnorm(n    = length(simpost_ccg),
                                                        mean = x["post_m"],
                                                        sd   = x["post_s"])})

        #-----------------------------------------------------------------------
        # Matrix with samples from posteriors

        simpost <- cbind(simpost_ccg,
                         simpost_trt)

        #-----------------------------------------------------------------------
        # Output object for weights

        simdat$w_post[nrow(simdat)]

        weights <- data.frame(w_prior = c(1 - robust, robust),
                              w_post  = c(1 - simdat$w_post[nrow(simdat)],
                                          simdat$w_post[nrow(simdat)]))

        #-----------------------------------------------------------------------
        # Prior parameters + estimates output

        m_est <- simdat$m
        s_est <- simdat$s

        estimates_out <- simdat

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
        # Effective posterior sample size

        eff_n_post <- try(round(ess(posterior, method = "elir")), silent = TRUE)

        if(inherits(eff_n_post, "try-error")){

                eff_n_post <- round(ess(posterior, method = "moment"))

                warning("Estimation of ESS via elir-method does not work.\n ESS was estimated based on method of moments")
        }

        eff_n_post <- min(max(0, eff_n_post),
                          sum(histdat$n))

        #-----------------------------------------------------------------------
        # Shrinkage

        postmean <- summary(posterior)["mean"]
        ccg_mean <- newdat$mean[1]

        shrink_abs <- unname(abs(ccg_mean - postmean))
        shrink_rel <- if(ccg_mean == 0) NA_real_ else shrink_abs / abs(ccg_mean)
        shrinkage  <- c(abs = shrink_abs, rel = shrink_rel)

        #-----------------------------------------------------------------------
        ## ratio of means (mean_1 / mean_2)

        if(contr == "mean_ratio"){

                if(any(newdat$mean <= 0)){
                        warning("Ratio of means is only meaningful for strictly positive means.")
                }

                # Get contrast matrix
                nni <- newdat$n
                names(nni) <- newdat$group

                cmatr <- mratios::contrMatRatio(n    = nni,
                                                type = type,
                                                base = base)

                #---------------------------------------------------------------
                # Compute ratios

                ratio2control <- BSagri::CCRatio.default(x = simpost, cmat = cmatr)
                dpostratio    <- as.data.frame(ratio2control$chains)

                #---------------------------------------------------------------
                # Simultaneous confidence sets from empirical joint distribution

                SCIR <- BSagri::SCSnp(x = ratio2control, conf.level = 1 - alpha,
                                      alternative = alt_scs)

                dscir <- data.frame(
                        row.names = seq_along(SCIR$estimate),
                        comp = names(SCIR$estimate),
                        estimate  = unname(signif(SCIR$estimate, 3)),
                        lower.ci   = unname(signif(SCIR$conf.int[,1], 3)),
                        upper.ci   = unname(signif(SCIR$conf.int[,2], 3)))
        }

        #-----------------------------------------------------------------------
        ## difference of means (mean_1 - mean_2)

        if(contr == "mean_diff"){

                # Get contrast matrix
                nni <- newdat$n
                names(nni) <- newdat$group

                cmatr <- multcomp::contrMat(n    = nni,
                                            type = type,
                                            base = base)

                #---------------------------------------------------------------
                # Compute differences

                diff2control <- BSagri::CCDiff.default(x = simpost, cmat = cmatr)
                dpostdiff    <- as.data.frame(diff2control$chains)

                #---------------------------------------------------------------
                # Simultaneous confidence sets from empirical joint distribution

                SCIR <- BSagri::SCSnp(x = diff2control, conf.level = 1 - alpha,
                                      alternative = alt_scs)

                
                dscir <- data.frame(
                        row.names = seq_along(SCIR$estimate),
                        comp = names(SCIR$estimate),
                        estimate  = unname(signif(SCIR$estimate, 3)),
                        lower.ci   = unname(signif(SCIR$conf.int[,1], 3)),
                        upper.ci   = unname(signif(SCIR$conf.int[,2], 3)))
        }

        #-----------------------------------------------------------------------
        # Output

        out_list <- list("comp"       = dscir,
                         "contr"      = contr,
                         "histdat"    = histdat,
                         "newdat"     = newdat,
                         "distr"      = "normal",
                         "prior_par"  = c("m" = m_est, "s" = s_est),
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
                            class = c("SCI", "MultiMAP_norm", "MultiMAP"))

        return(out_s3)
}


#-------------------------------------------------------------------------------
# Internal: run the borrowing pipeline for ONE arm (normal endpoint).
# 'hd' is that arm's historical data (mean, sd, n); (m_cur, s_cur, n_cur) the
# arm's current-trial summary. Returns the per-arm borrowing bundle, including
# n_mcmc posterior draws in $samples. Mirrors the single-arm control block.
#-------------------------------------------------------------------------------
.borrow_norm_arm <- function(hd, m_cur, s_cur, n_cur,
                             robust, sigma, tau_prior, tau_dist, beta_prior,
                             limits_pi, alpha, n_mcmc){

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

        # Prediction interval for a future group of this arm's size
        pred_dist <- RBesT::preddist(map, n = n_cur, sigma = sigma_ref)
        pred_int  <- .pred_int_from(pred_dist, limits_pi, alpha)

        # Robustify + posterior update with the arm's current data
        map_robust <- robustify(map, weight = robust,
                                mean = summary(map)["mean"], sigma = sigma_ref)
        posterior  <- postmix(map_robust, m = m_cur, se = s_cur / sqrt(n_cur))
        samples    <- RBesT::rmix(posterior, n_mcmc)

        cap        <- sum(hd$n)
        postmean   <- unname(summary(posterior)["mean"])
        shrink_abs <- abs(m_cur - postmean)
        shrink_rel <- if(m_cur == 0) NA_real_ else shrink_abs / abs(m_cur)
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
             tau_prior  = .tau_prior_out(tau_dist, tau_p),
             beta_prior = beta_p,
             samples    = samples)
}

#-------------------------------------------------------------------------------
# Internal: simultaneous credible intervals of the treatment-vs-control contrasts
# from the posterior-sample matrix (normal endpoint). Columns of 'simpost' are in
# newdat row order. Returns list(dscir, cmatr).
#-------------------------------------------------------------------------------
.sci_norm <- function(simpost, nni, contr, type, base, alpha, alt_scs){
        if(contr == "mean_ratio"){
                cmatr <- mratios::contrMatRatio(n = nni, type = type, base = base)
                cc    <- BSagri::CCRatio.default(x = simpost, cmat = cmatr)
        } else {
                cmatr <- multcomp::contrMat(n = nni, type = type, base = base)
                cc    <- BSagri::CCDiff.default(x = simpost, cmat = cmatr)
        }
        SCIR  <- BSagri::SCSnp(x = cc, conf.level = 1 - alpha, alternative = alt_scs)
        dscir <- data.frame(row.names = seq_along(SCIR$estimate),
                            comp     = names(SCIR$estimate),
                            estimate = unname(signif(SCIR$estimate, 3)),
                            lower.ci = unname(signif(SCIR$conf.int[, 1], 3)),
                            upper.ci = unname(signif(SCIR$conf.int[, 2], 3)))
        list(dscir = dscir, cmatr = cmatr)
}

#-------------------------------------------------------------------------------
# Internal: multi-arm SCI_norm(). Each arm named in 'histdat' is borrowed (its
# own MAP prior + posterior); every other arm uses the non-borrowing flat-prior
# posterior. Contrasts and SCIs are formed exactly as in the single-arm path.
#-------------------------------------------------------------------------------
.sci_norm_multiarm <- function(histdat, newdat, contr, type, base,
                               robust, sigma, tau_prior, tau_dist, beta_prior,
                               limits, limits_pi, alpha, n_mcmc, seed){

        arms        <- as.character(newdat$group)
        borrow_arms <- names(histdat)

        alt_scs <- switch(limits, "lower" = "greater", "upper" = "less", "two.sided")

        # Per-arm posterior samples (newdat row order); borrowing bundles by arm.
        cols   <- vector("list", length(arms))
        borrow <- vector("list", 0)
        for(i in seq_along(arms)){
                a <- arms[i]
                if(a %in% borrow_arms){
                        b <- .borrow_norm_arm(histdat[[a]],
                                              newdat$mean[i], newdat$sd[i], newdat$n[i],
                                              robust, sigma, tau_prior, tau_dist, beta_prior,
                                              limits_pi, alpha, n_mcmc)
                        cols[[i]]   <- b$samples
                        borrow[[a]] <- b
                } else {
                        # Non-borrowing arm: improper prior on the mean ->
                        # mean ~ N(ybar, sd / sqrt(n)).
                        cols[[i]] <- rnorm(n_mcmc, mean = newdat$mean[i],
                                           sd = newdat$sd[i] / sqrt(newdat$n[i]))
                }
        }
        simpost <- do.call(cbind, cols)
        colnames(simpost) <- arms

        if(contr == "mean_ratio" && any(newdat$mean <= 0))
                warning("Ratio of means is only meaningful for strictly positive means.")

        nni        <- newdat$n; names(nni) <- newdat$group
        sci        <- .sci_norm(simpost, nni, contr, type, base, alpha, alt_scs)

        # Keep the borrowing bundles in newdat order for tidy per-arm output.
        borrow <- borrow[intersect(arms, names(borrow))]
        pick   <- function(nm) lapply(borrow, `[[`, nm)
        pickv  <- function(nm) vapply(borrow, `[[`, numeric(1), nm)

        out_list <- list("comp"        = sci$dscir,
                         "contr"       = contr,
                         "histdat"     = histdat,
                         "newdat"      = newdat,
                         "distr"       = "normal",
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

        structure(out_list, class = c("SCI", "MultiMAP_norm", "MultiMAP"))
}


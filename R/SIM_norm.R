#' Monte-Carlo operating characteristics for the normal MultiMAP workflow
#'
#' Runs a complete Monte-Carlo study for one parameter scenario and returns a
#' one-row data frame of operating characteristics, so rows for many scenarios
#' can simply be \code{rbind()}ed (e.g. over an \code{expand.grid()}). For each
#' of the \code{n_sim} replicates it samples \code{H} historical control
#' studies and one current trial from \code{\link{h1_norm}} (each study's
#' within-study residual SD drawn from a gamma), fits the borrowing
#' analysis \code{\link{SCI_norm}}, and records the per-contrast reject
#' decisions, the prediction-interval decision and the borrowing diagnostics;
#' these are then averaged over replicates.
#'
#' @param n_sim number of Monte-Carlo replicates (simulated trials).
#' @param H number of historical control studies.
#' @param n number of observations per arm (control and each treatment).
#' @param mu grand mean of the data-generating process (shared by historical and
#'   current controls -- i.e. no prior-data conflict).
#' @param sd_h between-study (random effect) standard deviation.
#' @param sd_hi_mean mean of the gamma distribution from which each study's
#'   within-study residual SD is drawn.
#' @param sd_hi_cv coefficient of variation of that gamma distribution
#'   (\code{sd_hi_cv = 0} gives a fixed within-study SD equal to
#'   \code{sd_hi_mean}). The gamma has shape \code{1 / sd_hi_cv^2} and rate
#'   \code{shape / sd_hi_mean}.
#' @param shift true treatment-vs-control effect(s); \code{length(shift)} fixes
#'   the number of treatment groups. All-zero -> null (FWER); non-zero ->
#'   alternative (power).
#' @param mu_cur grand mean of the CURRENT trial. \code{NULL} (default) =
#'   \code{mu}, i.e. no prior-data conflict. A value != \code{mu} shifts the
#'   whole current trial (concurrent control AND every treatment arm) by
#'   \code{mu_cur - mu} relative to the historical controls, inducing a
#'   prior-data conflict while leaving the treatment-vs-control effects unchanged.
#' @param prior_est analysis for each replicate, passed to \code{\link{SCI_norm}}:
#'   \code{"MAP"} (default, dynamic borrowing) or \code{"ignore_HCD"} (the
#'   non-borrowing frequentist competitor). Under \code{"ignore_HCD"} the
#'   borrowing-specific output columns (\code{pi_t1e}, \code{pi_power},
#'   \code{shrink_*}, \code{eff_n_*}, \code{post_w_*}, ...) are \code{NA}, so the
#'   two runs can be \code{rbind()}ed.
#' @param type multiple-comparison type passed through to \code{\link{SCI_norm}}:
#'   either \code{"Dunnett"} (many-to-one, each treatment vs. control; the
#'   default) or \code{"Tukey"} (all pairwise comparisons). Under \code{"Tukey"}
#'   the per-contrast \code{reject_*} output columns are keyed by the full
#'   pairwise contrast labels.
#' @param contr \code{"mean_diff"} (differences) or \code{"mean_ratio"} (ratios
#'   of means).
#' @param limits \code{"two.sided"} (default), \code{"lower"} or
#'   \code{"upper"} (drives the reject decision; passed through to
#'   \code{\link{SCI_norm}}).
#' @param base index of the control group in the contrast matrix.
#' @param robust weight of the weakly-informative (robust) prior component.
#' @param sigma reference scale for the RBesT Gaussian machinery; \code{NULL}
#'   estimates it from the historical controls.
#' @param tau_prior,beta_prior optional overrides for the \code{gMAP} between-study
#'   heterogeneity and intercept priors, passed to \code{\link{SCI_norm}}.
#'   \code{NULL} (the default) uses the built-in heuristics. Increasing the
#'   \code{tau_prior} scale widens the MAP prior and can improve small-\code{H}
#'   prediction-interval calibration.
#' @param alpha 1 - simultaneous confidence level.
#' @param n_mcmc number of Monte-Carlo draws inside each borrowing fit.
#' @param parallel \code{TRUE} (default) fits the replicates on a PSOCK cluster
#'   (works on Windows, macOS and Linux); \code{FALSE} runs them serially.
#' @param n_cores number of worker processes when \code{parallel = TRUE}.
#'   \code{NULL} (default) uses \code{parallel::detectCores() - 1}, capped at
#'   \code{n_sim}.
#' @param seed optional integer. When supplied, one independent L'Ecuyer-CMRG RNG
#'   stream is assigned per replicate, so the study is reproducible AND identical
#'   whether run serially or in parallel and regardless of the number of workers.
#'   \code{NULL} gives a non-reproducible run.
#'
#' @return A one-row data frame of operating characteristics, including:
#'   \describe{
#'     \item{\code{fwer}}{family-wise error rate = P(>= 1 contrast rejected) under
#'       the global null (all shift == 0); \code{NA} otherwise.}
#'     \item{\code{power}}{any-pairs power = P(>= 1 contrast rejected) when at
#'       least one shift != 0; \code{NA} under the global null. The per-contrast
#'       decision follows \code{limits}: \code{"lower"} lower limit only,
#'       \code{"upper"} upper limit only, \code{"two.sided"} either limit.}
#'     \item{\code{pi_t1e}}{type-1 error of the MAP-based prediction interval (the
#'       rate the concurrent control mean falls outside it), reported only under
#'       NO prior-data conflict (\code{mu_cur == mu}); \code{NA} otherwise.}
#'     \item{\code{pi_power}}{the same rate reported only WITH conflict
#'       (\code{mu_cur != mu}), i.e. the PI's power to detect the whole-trial
#'       shift; \code{NA} under no conflict. (pi_t1e / pi_power mirror the
#'       fwer / power split.)}
#'     \item{\code{pi_reject_lower}, \code{pi_reject_upper}}{PI rejection rate
#'       split by border (always reported).}
#'     \item{\code{shrink_abs}, \code{shrink_rel}}{mean absolute / relative
#'       shrinkage of the concurrent control mean towards the borrowed prior.}
#'     \item{\code{width_mean}, \code{eff_n_mean}, \code{eff_n_post_mean},
#'       \code{post_w_mean}, \code{post_w_median}}{further borrowing diagnostics
#'       averaged over replicates.}
#'     \item{\code{reject_*}, \code{reject_lower_*}, \code{reject_upper_*}}{
#'       per-contrast rejection rates.}
#'   }
#'
#' @seealso \code{\link{SCI_norm}}, \code{\link{h1_norm}}
#'
#' @examples
#' \donttest{
#' # Null scenario -> FWER and PI type-1 error
#' SIM_norm(n_sim = 100, H = 5, n = 30, mu = 10, sd_h = 2,
#'          sd_hi_mean = 2, sd_hi_cv = 0.5, shift = c(0, 0),
#'          seed = 1, parallel = TRUE, n_cores = 2)
#'
#' # Prior-data conflict -> PI power to detect the shift
#' SIM_norm(n_sim = 100, H = 5, n = 30, mu = 10, sd_h = 2,
#'          sd_hi_mean = 2, sd_hi_cv = 0.5, shift = c(0, 1.5),
#'          mu_cur = 13, seed = 1, parallel = TRUE, n_cores = 2)
#'
#' # Non-borrowing frequentist competitor (unequal-variance Dunnett) on the same
#' # simulated trials -- rbind() with the borrowing run to compare. The same
#' # seed keeps the current trials identical, so the comparison is paired.
#' SIM_norm(n_sim = 100, H = 5, n = 30, mu = 10, sd_h = 2,
#'          sd_hi_mean = 2, sd_hi_cv = 0.5, shift = c(0, 0),
#'          prior_est = "ignore_HCD", seed = 1, parallel = TRUE, n_cores = 2)
#' }
#' @export
SIM_norm <- function(n_sim,
                         H,
                         n,
                         mu,
                         sd_h,
                         sd_hi_mean,
                         sd_hi_cv,
                         shift       = c(0, 0),
                         mu_cur      = NULL,
                         prior_est   = "MAP",
                         type        = "Dunnett",
                         contr       = "mean_diff",
                         limits = "two.sided",
                         base        = 1,
                         robust      = 0.2,
                         sigma       = NULL,
                         tau_prior   = NULL,
                         beta_prior  = NULL,
                         alpha       = 0.05,
                         n_mcmc      = 50000,
                         parallel    = TRUE,
                         n_cores     = NULL,
                         seed        = NULL){

        # Number of treatment groups is defined by length(shift)
        m           <- length(shift)
        null_val    <- .null_value(contr)
        limits <- match.arg(limits, c("two.sided", "lower", "upper"))
        # Multiple-comparison type: many-to-one (Dunnett) or all-pairwise (Tukey).
        type        <- match.arg(type, c("Dunnett", "Tukey"))
        if(!prior_est %in% c("MAP", "ignore_HCD"))
                stop("'prior_est' must be 'MAP' or 'ignore_HCD'.")
        if(prior_est == "ignore_HCD" && contr != "mean_diff")
                stop("prior_est = 'ignore_HCD' supports only contr = 'mean_diff' for the normal endpoint.")

        # Current-trial grand mean. NULL -> mu (no prior-data conflict). A value
        # != mu shifts the WHOLE current trial (concurrent control AND every
        # treatment arm) by mu_cur - mu relative to the historical controls, i.e.
        # induces a prior-data conflict while leaving the treatment-vs-control
        # effects ('shift') unchanged.
        if(is.null(mu_cur)) mu_cur <- mu

        # Force all arguments before the parallel dispatch, so run_one()'s closure
        # captures their VALUES rather than promises bound to the caller's
        # environment. Without this, PSOCK workers cannot resolve an argument
        # passed as an unevaluated variable expression (e.g. H = H).
        for(.nm in names(formals())) force(get(.nm, envir = environment(), inherits = FALSE))

        # Silence RBesT / rstan chatter; genuine errors still propagate.
        quiet <- function(expr) suppressWarnings(suppressMessages(expr))

        # Draw k within-study residual SDs from a gamma with mean sd_hi_mean and
        # coefficient of variation sd_hi_cv (cv = 0 -> a fixed value).
        draw_sd_hi <- function(k){
                if(sd_hi_cv <= 0) return(rep(sd_hi_mean, k))
                shp <- 1 / sd_hi_cv^2
                rgamma(k, shape = shp, rate = shp / sd_hi_mean)
        }

        #-----------------------------------------------------------------------
        # Per-replicate operating characteristics from one fitted SCI_norm() object

        op_char <- function(fit){

                ccg_mean  <- fit$newdat$mean[1]
                post_mean <- sum(fit$posterior["w", ] * fit$posterior["m", ])

                # Per-contrast reject decisions, resolved by credible limit
                reject_lower <- fit$comp$lower.ci > null_val
                reject_upper <- fit$comp$upper.ci < null_val
                width        <- fit$comp$upper.ci - fit$comp$lower.ci

                # Prediction-interval decision for the concurrent control mean
                reject_pl <- unname(fit$pred_int["lower"] > ccg_mean)
                reject_pu <- unname(fit$pred_int["upper"] < ccg_mean)

                # Control shrinkage (absolute + relative)
                shrink_abs <- abs(ccg_mean - post_mean)
                shrink_rel <- shrink_abs / abs(ccg_mean)

                list(reject_lower = reject_lower,
                     reject_upper = reject_upper,
                     width        = width,
                     reject_pl    = reject_pl,
                     reject_pu    = reject_pu,
                     reject_pi    = isTRUE(reject_pl) || isTRUE(reject_pu),
                     shrink_abs   = shrink_abs,
                     shrink_rel   = shrink_rel,
                     eff_n        = fit$eff_n,
                     eff_n_post   = fit$eff_n_post,
                     w_post       = fit$weights$w_post[1],
                     trt_names    = as.character(fit$newdat$group[-1]),
                     comp_names   = as.character(fit$comp$comp))
        }

        #-----------------------------------------------------------------------
        # Reproducibility. When 'seed' is supplied, one independent L'Ecuyer-CMRG
        # RNG stream is assigned per replicate; setting it inside run_one() makes
        # the results identical whether run serially or in parallel and for any
        # number of workers. The caller's RNG kind/state is restored on exit.

        if(!is.null(seed)){
                if(exists(".Random.seed", envir = .GlobalEnv)){
                        .old_seed <- get(".Random.seed", envir = .GlobalEnv)
                        on.exit(assign(".Random.seed", .old_seed, envir = .GlobalEnv), add = TRUE)
                }
                .old_kind <- RNGkind()
                on.exit(suppressWarnings(do.call(RNGkind, as.list(.old_kind))), add = TRUE)
        }

        make_streams <- function(seed, n){
                if(is.null(seed)) return(vector("list", n))   # NULL -> no seeding
                RNGkind("L'Ecuyer-CMRG")
                set.seed(seed)
                s      <- vector("list", n)
                s[[1]] <- get(".Random.seed", envir = .GlobalEnv)
                for(i in seq_len(n)[-1]) s[[i]] <- parallel::nextRNGStream(s[[i - 1]])
                s
        }
        streams <- make_streams(seed, n_sim)

        #-----------------------------------------------------------------------
        # One Monte-Carlo replicate: sample HCD + current trial, fit SCI_norm(),
        # reduce to operating characteristics. Returns NULL if the fit fails.

        run_one <- function(stream){

                if(!is.null(stream)){
                        RNGkind("L'Ecuyer-CMRG")
                        assign(".Random.seed", stream, envir = .GlobalEnv)
                }

                # Historical control studies (each study's within-study SD ~ gamma)
                sd_hi_hist <- draw_sd_hi(H)
                hist_raw   <- h1_norm(mu = mu, H = H, n = n, sd_h = sd_h,
                                 sd_hi = sd_hi_hist, shift = NULL)
                histdat    <- do.call(rbind, lapply(split(hist_raw$y_hi, hist_raw$a),
                                        function(y) data.frame(mean = mean(y),
                                                               sd   = sd(y),
                                                               n    = length(y))))

                # Current trial. mu_cur shifts the whole trial vs. the HCD (conflict).
                sd_hi_cur <- draw_sd_hi(1)
                cur_raw   <- h1_norm(mu = mu_cur, H = 1, n = n, sd_h = sd_h,
                                sd_hi = sd_hi_cur, shift = shift)
                newdat    <- do.call(rbind, lapply(levels(cur_raw$arm), function(g){
                                        y <- cur_raw$y_hi[cur_raw$arm == g]
                                        data.frame(arm  = g, mean = mean(y),
                                                   sd = sd(y), n = length(y))
                                }))

                # Non-borrowing frequentist competitor: analyse the current trial
                # alone with an unequal-variance Dunnett-type multiple comparison
                # (no HCD used, no MAP). Borrowing-specific quantities are NA.
                if(prior_est == "ignore_HCD"){
                        fc <- try(quiet(.freq_competitor(newdat, distr = "normal", type = type,
                                                         base = base, alpha = alpha,
                                                         alternative = limits)),
                                  silent = TRUE)
                        if(inherits(fc, "try-error")) return(NULL)
                        cmp <- fc$comp
                        return(list(reject_lower = cmp$lower.ci > null_val,
                                    reject_upper = cmp$upper.ci < null_val,
                                    width        = cmp$upper.ci - cmp$lower.ci,
                                    reject_pl    = NA, reject_pu = NA, reject_pi = NA,
                                    shrink_abs   = NA_real_, shrink_rel = NA_real_,
                                    eff_n        = NA_real_, eff_n_post = NA_real_, w_post = NA_real_,
                                    trt_names    = as.character(newdat[[1]][-1]),
                                    comp_names   = as.character(cmp$comp)))
                }

                fit <- try(quiet(SCI_norm(histdat = histdat, newdat = newdat,
                                          prior_est = prior_est, type = type, contr = contr,
                                          limits = limits, base = base,
                                          robust = robust, sigma = sigma,
                                          tau_prior = tau_prior, beta_prior = beta_prior,
                                          alpha = alpha, n_mcmc = n_mcmc)),
                           silent = TRUE)
                if(inherits(fit, "try-error")) return(NULL)
                op_char(fit)
        }

        #-----------------------------------------------------------------------
        # Fit all replicates, serially or on a PSOCK cluster. PSOCK works on
        # Windows, macOS and Linux; the per-replicate L'Ecuyer streams make the
        # serial and parallel paths identical for any number of workers.

        if(is.null(n_cores)){
                nc      <- parallel::detectCores()
                n_cores <- if(is.na(nc)) 1L else max(1L, nc - 1L)
        }
        n_cores <- min(as.integer(n_cores), n_sim)

        if(isTRUE(parallel) && n_cores > 1L){

                cl <- parallel::makeCluster(n_cores)          # PSOCK (all OSes)
                on.exit(parallel::stopCluster(cl), add = TRUE)

                # analysis packages on every worker (SCI_norm() needs these)
                parallel::clusterEvalQ(cl, suppressWarnings(suppressPackageStartupMessages({
                        library(RBesT); library(multcomp)
                        library(mratios); library(BSagri)
                })))
                # mirror the caller's RBesT MCMC options; run Stan chains
                # sequentially inside each worker to avoid oversubscription
                rbest_opts <- options()[grep("^RBesT\\.", names(options()))]
                parallel::clusterCall(cl, function(o){ options(o); options(mc.cores = 1L) },
                                      rbest_opts)
                # ship the analysis functions (all scalars/params travel with
                # run_one()'s closure environment)
                parallel::clusterExport(cl, varlist = c("SCI_norm", "h1_norm", ".freq_competitor"),
                                        envir = environment())

                oc_list <- parallel::parLapply(cl, streams, run_one)
        } else {
                oc_list <- lapply(streams, run_one)
        }

        # Drop replicates whose borrowing fit failed (SCI_norm() wrapped in try())
        oc_list   <- Filter(Negate(is.null), oc_list)
        n_sim_eff <- length(oc_list)
        if(n_sim_eff < n_sim*0.75) warning("More than 25% of the replicates failed for this scenario.")
        if(n_sim_eff == 0) stop("All replicates failed for this scenario.")
        
        #=======================================================================
        # Aggregate the per-replicate results into the operating characteristics
        #=======================================================================
        
        # Contrast labels + count. For Dunnett there is one contrast per treatment
        # group, so the short treatment-group labels are kept (unchanged output);
        # for Tukey (all pairwise comparisons) there are choose(m + 1, 2) contrasts,
        # so the full contrast labels from the fit are used.
        trt_names   <- oc_list[[1]]$trt_names
        comp_names  <- oc_list[[1]]$comp_names
        n_contr     <- length(comp_names)
        contr_names <- if(identical(type, "Dunnett")) trt_names else comp_names

        rej_lower_mat <- t(vapply(oc_list, function(r) as.logical(r$reject_lower), logical(n_contr)))
        rej_upper_mat <- t(vapply(oc_list, function(r) as.logical(r$reject_upper), logical(n_contr)))
        
        reject_mat <- switch(limits,
                             "lower" = rej_lower_mat,
                             "upper" = rej_upper_mat,
                             rej_lower_mat | rej_upper_mat)   # two.sided
        
        any_row    <- function(mat) if(ncol(mat) > 0) apply(mat, 1, any) else logical(nrow(mat))
        any_reject <- mean(any_row(reject_mat))
        
        fwer  <- if(all(shift == 0)) any_reject else NA_real_
        power <- if(any(shift != 0)) any_reject else NA_real_
        
        fw_reject_lower <- mean(any_row(rej_lower_mat))
        fw_reject_upper <- mean(any_row(rej_upper_mat))
        
        reject       <- colMeans(reject_mat);    names(reject)       <- contr_names
        reject_lower <- colMeans(rej_lower_mat); names(reject_lower) <- contr_names
        reject_upper <- colMeans(rej_upper_mat); names(reject_upper) <- contr_names
        
        # Prediction-interval type-1 error (and per border)
        # Rate at which the concurrent control mean falls outside the MAP-based
        # prediction interval. Under NO prior-data conflict (mu_cur == mu) this is
        # the PI's type-1 error; WITH conflict (mu_cur != mu) it is the PI's power
        # to detect the shift -- mirroring the fwer / power split for the SCI.
        pi_reject       <- mean(vapply(oc_list, function(r) r$reject_pi, logical(1)))
        pi_reject_lower <- mean(vapply(oc_list, function(r) r$reject_pl, logical(1)))
        pi_reject_upper <- mean(vapply(oc_list, function(r) r$reject_pu, logical(1)))

        pi_conflict <- !isTRUE(all.equal(mu_cur, mu))          # whole-trial shift vs HCD
        pi_t1e   <- if(!pi_conflict) pi_reject else NA_real_   # no conflict  -> type-1 error
        pi_power <- if( pi_conflict) pi_reject else NA_real_   # conflict     -> power to detect it
        
        # Shrinkage + further diagnostics
        shrink_abs      <- mean(vapply(oc_list, function(r) r$shrink_abs, numeric(1)))
        shrink_rel      <- mean(vapply(oc_list, function(r) r$shrink_rel, numeric(1)))
        width_mean      <- mean(vapply(oc_list, function(r) mean(r$width), numeric(1)))
        eff_n_mean      <- mean(vapply(oc_list, function(r) r$eff_n,      numeric(1)))
        eff_n_post_mean <- mean(vapply(oc_list, function(r) r$eff_n_post, numeric(1)))
        post_w          <- vapply(oc_list, function(r) r$w_post, numeric(1))
        post_w_mean     <- mean(post_w)
        post_w_median   <- median(post_w)
        
        #-----------------------------------------------------------------------
        # Assemble the one-row output. tau_prior / beta_prior are recorded as a
        # compact string ("default" when NULL, else the comma-joined values).

        tau_prior_str  <- if(is.null(tau_prior))  "default" else paste(format(tau_prior),  collapse = ",")
        beta_prior_str <- if(is.null(beta_prior)) "default" else paste(format(beta_prior), collapse = ",")

        out <- data.frame(n_sim            = n_sim,
                          n_sim_eff        = n_sim_eff,
                          H               = H,
                          n               = n,
                          mu              = mu,
                          mu_cur          = mu_cur,
                          sd_h            = sd_h,
                          sd_hi_mean      = sd_hi_mean,
                          sd_hi_cv        = sd_hi_cv,
                          m               = m,
                          shift           = paste(shift, collapse = ";"),
                          prior_est       = prior_est,
                          type            = type,
                          contr           = contr,
                          limits     = limits,
                          base            = base,
                          robust          = robust,
                          sigma           = if(is.null(sigma)) NA_real_ else sigma,
                          tau_prior       = tau_prior_str,
                          beta_prior      = beta_prior_str,
                          alpha           = alpha,
                          n_mcmc          = n_mcmc,
                          seed            = if(is.null(seed)) NA_integer_ else seed,
                          
                          fwer            = fwer,    # global null only, else NA
                          power           = power,   # alternative only, else NA
                          
                          fw_reject_lower = fw_reject_lower,
                          fw_reject_upper = fw_reject_upper,
                          
                          pi_t1e          = pi_t1e,    # no-conflict scenarios only, else NA
                          pi_power        = pi_power,  # conflict scenarios only, else NA
                          pi_reject_lower = pi_reject_lower,
                          pi_reject_upper = pi_reject_upper,
                          
                          shrink_abs      = shrink_abs,
                          shrink_rel      = shrink_rel,
                          
                          width_mean      = width_mean,
                          eff_n_mean      = eff_n_mean,
                          eff_n_post_mean = eff_n_post_mean,
                          post_w_mean     = post_w_mean,
                          post_w_median   = post_w_median,
                          stringsAsFactors = FALSE)
        
        # Per-contrast rejection rates, one column block per contrast
        for(g in contr_names){
                out[[paste0("reject_",       g)]] <- reject[[g]]
                out[[paste0("reject_lower_", g)]] <- reject_lower[[g]]
                out[[paste0("reject_upper_", g)]] <- reject_upper[[g]]
        }
        
        out
}

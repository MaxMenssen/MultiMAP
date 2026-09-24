#' Monte-Carlo operating characteristics for the Poisson MultiMAP workflow
#'
#' Runs a complete Monte-Carlo study for one parameter scenario and returns a
#' one-row data frame of operating characteristics, so rows for many scenarios
#' can simply be \code{rbind()}ed (e.g. over an \code{expand.grid()}). The
#' Poisson analogue of \code{\link{SIM_norm}} and \code{\link{SIM_binom}}. For
#' each of the \code{n_sim} replicates it samples \code{H} historical control
#' studies and one current trial from \code{\link{h1_pois}} (a gamma-Poisson
#' random-effects model), fits the borrowing analysis \code{\link{SCI_pois}},
#' and records the per-contrast reject decisions, the prediction-interval
#' decision and the borrowing diagnostics; these are then averaged over
#' replicates.
#'
#' @section Between-study variation:
#' The gamma-Poisson between-study variation is set with exactly one of
#' \code{shape}, \code{sd_h} or \code{cv_h}. \code{shape} is the gamma shape of
#' the between-study rate distribution (Tarone's \code{p}); the between-study
#' variance of the true cluster rates is \eqn{\lambda^2 / \mathrm{shape}}, so a
#' between-study standard deviation \code{sd_h} or coefficient of variation
#' \code{cv_h} is converted to \code{shape} via
#' \eqn{\mathrm{shape} = \lambda^2 / \mathrm{sd\_h}^2} and
#' \eqn{\mathrm{shape} = 1 / \mathrm{cv\_h}^2} (using the historical
#' \code{lambda}). Larger \code{shape} means less heterogeneity.
#'
#' @param n_sim number of Monte-Carlo replicates (simulated trials).
#' @param H number of historical control studies.
#' @param offset exposure (number of experimental units) per arm.
#' @param lambda grand rate of the data-generating process (shared by historical
#'   and current controls -- i.e. no prior-data conflict).
#' @param shape gamma shape of the between-study rate distribution (between-study
#'   variation). Supply exactly one of \code{shape}, \code{sd_h}, \code{cv_h}.
#' @param sd_h between-study standard deviation of the true cluster rates
#'   (converted to \code{shape}; see the between-study variation section).
#' @param cv_h between-study coefficient of variation of the true cluster rates
#'   (converted to \code{shape}).
#' @param shift true treatment-vs-control effect(s) on the scale given by
#'   \code{shift_scale}; \code{length(shift)} fixes the number of treatment
#'   groups. All-zero -> null (FWER); non-zero -> alternative (power).
#' @param shift_scale scale on which \code{shift} acts: \code{"rate"} (default,
#'   rate difference) or \code{"log"} (log rate ratio). See \code{\link{h1_pois}}.
#' @param lambda_cur grand rate of the CURRENT trial. \code{NULL} (default) =
#'   \code{lambda}, i.e. no prior-data conflict. A value != \code{lambda} shifts
#'   the whole current trial (concurrent control AND every treatment arm)
#'   relative to the historical controls, inducing a prior-data conflict while
#'   leaving the treatment-vs-control effects unchanged.
#' @param prior_est analysis for each replicate: \code{"MAP"} (default) or
#'   \code{"EB"} borrowing (passed to \code{\link{SCI_pois}}), or
#'   \code{"ignore_HCD"} -- a non-borrowing frequentist competitor fitted here
#'   from the current trial alone (Poisson GLM with a log offset + unequal-
#'   variance Dunnett-type multiple comparison on rate ratios; requires
#'   \code{contr = "rate_ratio"}). Under \code{"ignore_HCD"} the borrowing-specific
#'   output columns (\code{pi_t1e}, \code{pi_power}, \code{shrink_*},
#'   \code{eff_n_*}, \code{post_w_*}, ...) are \code{NA}, so the runs
#'   \code{rbind()} together.
#' @param type multiple-comparison type passed through to \code{\link{SCI_pois}}:
#'   either \code{"Dunnett"} (many-to-one, each treatment vs. control; the
#'   default) or \code{"Tukey"} (all pairwise comparisons). Under \code{"Tukey"}
#'   the per-contrast \code{reject_*} output columns are keyed by the full
#'   pairwise contrast labels.
#' @param contr \code{"rate_ratio"} (ratio of rates; default) or
#'   \code{"rate_diff"} (difference of rates).
#' @param limits \code{"two.sided"} (default), \code{"lower"} or
#'   \code{"upper"} (drives the reject decision; passed to \code{\link{SCI_pois}}).
#' @param base index of the control group in the contrast matrix.
#' @param robust weight of the weakly-informative (robust) prior component.
#' @param tau_prior,beta_prior optional overrides for the \code{gMAP} between-study
#'   heterogeneity and intercept priors (on the log-rate scale), passed to
#'   \code{\link{SCI_pois}} (MAP path only). \code{NULL} (the default) uses the
#'   built-in heuristics, whose reference scale is floored so the default
#'   \code{tau_prior} is a weakly-informative \code{HalfNormal} of scale at least
#'   1 on the log-rate scale (see \code{\link{SCI_pois}}).
#' @param alpha 1 - simultaneous confidence level.
#' @param n_mcmc number of Monte-Carlo draws inside each borrowing fit (passed to
#'   \code{\link{SCI_pois}}).
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
#' @return A one-row data frame of operating characteristics, including
#'   \code{fwer}, \code{power}, \code{pi_t1e} / \code{pi_power} (the
#'   prediction-interval type-1 error / power, split by whether
#'   \code{lambda_cur == lambda}), shrinkage and further borrowing diagnostics,
#'   and the per-contrast rejection rates. See \code{\link{SIM_binom}} for the
#'   column descriptions.
#'
#' @seealso \code{\link{SCI_pois}}, \code{\link{h1_pois}}, \code{\link{SIM_binom}}
#'
#' @examples
#' \donttest{
#' # Null scenario -> FWER and PI type-1 error
#' SIM_pois(n_sim = 100, H = 6, offset = 3, lambda = 8, shape = 12,
#'          shift = c(0, 0), seed = 1, parallel = TRUE, n_cores = 2)
#'
#' # Prior-data conflict -> PI power to detect the shift
#' SIM_pois(n_sim = 100, H = 6, offset = 3, lambda = 8, shape = 12,
#'          shift = c(0, 4), lambda_cur = 4, seed = 1,
#'          parallel = TRUE, n_cores = 2)
#'
#' # Non-borrowing frequentist competitor (Poisson GLM + unequal-variance Dunnett
#' # on rate ratios). rbind() with a borrowing run to compare; the shared seed
#' # keeps the current trials identical (paired comparison).
#' SIM_pois(n_sim = 100, H = 6, offset = 3, lambda = 8, shape = 12,
#'          shift = c(0, 0), prior_est = "ignore_HCD", seed = 1,
#'          parallel = TRUE, n_cores = 2)
#' }
#' @export
SIM_pois <- function(n_sim,
                     H,
                     offset,
                     lambda,
                     shape       = NULL,
                     sd_h        = NULL,
                     cv_h        = NULL,
                     shift       = c(0, 0),
                     shift_scale = c("rate", "log"),
                     lambda_cur  = NULL,
                     prior_est   = "MAP",
                     type        = "Dunnett",
                     contr       = "rate_ratio",
                     limits = "two.sided",
                     base        = 1,
                     robust      = 0.2,
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
        shift_scale <- match.arg(shift_scale)

        if(!is.numeric(lambda) || length(lambda) != 1 || lambda <= 0)
                stop("'lambda' must be a single positive number.")
        if(!prior_est %in% c("MAP", "EB", "ignore_HCD"))
                stop("'prior_est' must be 'MAP', 'EB' or 'ignore_HCD'.")
        if(prior_est == "ignore_HCD" && contr != "rate_ratio")
                stop("prior_est = 'ignore_HCD' supports only contr = 'rate_ratio' for the Poisson endpoint.")

        #-----------------------------------------------------------------------
        # Between-study variation: convert sd_h / cv_h to the canonical shape.
        # var_between(lambda_h) = lambda^2 / shape.

        n_supplied <- sum(!is.null(shape), !is.null(sd_h), !is.null(cv_h))
        if(n_supplied != 1)
                stop("Supply exactly one of 'shape', 'sd_h' or 'cv_h'.")

        if(!is.null(sd_h)) shape <- lambda^2 / sd_h^2
        if(!is.null(cv_h)) shape <- 1 / cv_h^2

        if(!is.numeric(shape) || length(shape) != 1 || shape <= 0)
                stop("The between-study variation is out of range: 'shape' must be positive.")

        # Between-study SD / CV of the true cluster rates (reported)
        sd_h_between <- lambda / sqrt(shape)
        cv_h_between   <- 1 / sqrt(shape)

        # Current-trial grand rate. NULL -> lambda (no prior-data conflict).
        if(is.null(lambda_cur)) lambda_cur <- lambda
        if(!is.numeric(lambda_cur) || length(lambda_cur) != 1 || lambda_cur <= 0)
                stop("'lambda_cur' must be a single positive number.")

        # Force all arguments before the parallel dispatch, so run_one()'s closure
        # captures their VALUES rather than promises bound to the caller's
        # environment. Without this, PSOCK workers cannot resolve an argument
        # passed as an unevaluated variable expression (e.g. H = H).
        for(.nm in names(formals())) force(get(.nm, envir = environment(), inherits = FALSE))

        # Silence RBesT / rstan chatter; genuine errors still propagate.
        quiet <- function(expr) suppressWarnings(suppressMessages(expr))

        #-----------------------------------------------------------------------
        # Per-replicate operating characteristics from one fitted SCI_pois() object

        op_char <- function(fit){

                ccg_mean  <- fit$newdat$mean[1]
                pw        <- fit$posterior["w", ]
                pa        <- fit$posterior["a", ]
                pb        <- fit$posterior["b", ]
                post_mean <- sum(pw * pa / pb)

                # Per-contrast reject decisions, resolved by credible limit
                reject_lower <- fit$comp$lower.ci > null_val
                reject_upper <- fit$comp$upper.ci < null_val
                width        <- fit$comp$upper.ci - fit$comp$lower.ci

                # Prediction-interval decision for the concurrent control rate
                reject_pl <- unname(fit$pred_int["lower"] > ccg_mean)
                reject_pu <- unname(fit$pred_int["upper"] < ccg_mean)

                # Control shrinkage (absolute + relative)
                shrink_abs <- abs(ccg_mean - post_mean)
                shrink_rel <- if(ccg_mean == 0) NA_real_ else shrink_abs / abs(ccg_mean)

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
        # Reproducibility (one independent L'Ecuyer-CMRG stream per replicate).

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
        # One Monte-Carlo replicate: sample HCD + current trial, fit SCI_pois(),
        # reduce to operating characteristics. Returns NULL if the fit fails.

        run_one <- function(stream){

                if(!is.null(stream)){
                        RNGkind("L'Ecuyer-CMRG")
                        assign(".Random.seed", stream, envir = .GlobalEnv)
                }

                # Historical control studies (control arm only, gamma-Poisson)
                hist_raw <- h1_pois(lambda = lambda, H = H, offset = offset,
                                    shape = shape, shift = NULL, shift_scale = shift_scale)
                histdat  <- data.frame(events = hist_raw$events,
                                       offset = hist_raw$offset)

                # Current trial (one study). lambda_cur shifts it vs. the HCD;
                # the treatment arms share the drawn cluster rate.
                cur_raw <- h1_pois(lambda = lambda_cur, H = 1, offset = offset,
                                   shape = shape, shift = shift, shift_scale = shift_scale)
                newdat  <- data.frame(group  = as.character(cur_raw$arm),
                                      events = cur_raw$events,
                                      offset = cur_raw$offset,
                                      stringsAsFactors = FALSE)

                # Non-borrowing frequentist competitor: Poisson GLM with a log
                # offset + unequal-variance Dunnett-type multiple comparison (rate
                # ratios), current trial only. Borrowing-specific quantities are NA.
                if(prior_est == "ignore_HCD"){
                        fc <- try(quiet(.freq_competitor(newdat, distr = "poisson", type = type,
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

                fit <- try(quiet(SCI_pois(histdat = histdat, newdat = newdat,
                                          prior_est = prior_est, type = type, contr = contr,
                                          limits = limits, base = base,
                                          robust = robust, tau_prior = tau_prior,
                                          beta_prior = beta_prior, alpha = alpha, n_mcmc = n_mcmc)),
                           silent = TRUE)
                if(inherits(fit, "try-error")) return(NULL)
                op_char(fit)
        }

        #-----------------------------------------------------------------------
        # Fit all replicates, serially or on a PSOCK cluster (all OSes).

        if(is.null(n_cores)){
                nc      <- parallel::detectCores()
                n_cores <- if(is.na(nc)) 1L else max(1L, nc - 1L)
        }
        n_cores <- min(as.integer(n_cores), n_sim)

        if(isTRUE(parallel) && n_cores > 1L){

                cl <- parallel::makeCluster(n_cores)          # PSOCK (all OSes)
                on.exit(parallel::stopCluster(cl), add = TRUE)

                parallel::clusterEvalQ(cl, suppressWarnings(suppressPackageStartupMessages({
                        library(RBesT); library(multcomp)
                        library(mratios); library(BSagri)
                })))
                rbest_opts <- options()[grep("^RBesT\\.", names(options()))]
                parallel::clusterCall(cl, function(o){ options(o); options(mc.cores = 1L) },
                                      rbest_opts)
                parallel::clusterExport(cl, varlist = c("SCI_pois", "h1_pois", ".freq_competitor"),
                                        envir = environment())

                oc_list <- parallel::parLapply(cl, streams, run_one)
        } else {
                oc_list <- lapply(streams, run_one)
        }

        # Drop replicates whose borrowing fit failed (SCI_pois() wrapped in try())
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

        # Prediction-interval type-1 error / power (split by conflict).
        pi_reject       <- mean(vapply(oc_list, function(r) r$reject_pi, logical(1)))
        pi_reject_lower <- mean(vapply(oc_list, function(r) r$reject_pl, logical(1)))
        pi_reject_upper <- mean(vapply(oc_list, function(r) r$reject_pu, logical(1)))

        pi_conflict <- !isTRUE(all.equal(lambda_cur, lambda))
        pi_t1e   <- if(!pi_conflict) pi_reject else NA_real_
        pi_power <- if( pi_conflict) pi_reject else NA_real_

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
                          offset          = offset,
                          lambda          = lambda,
                          lambda_cur      = lambda_cur,
                          shape           = shape,
                          sd_h_between    = sd_h_between,
                          cv_h_between      = cv_h_between,
                          m               = m,
                          shift           = paste(shift, collapse = ";"),
                          shift_scale     = shift_scale,
                          prior_est       = prior_est,
                          type            = type,
                          contr           = contr,
                          limits     = limits,
                          base            = base,
                          robust          = robust,
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

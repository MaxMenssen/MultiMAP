#' Monte-Carlo operating characteristics for the binomial MultiMAP workflow
#'
#' Runs a complete Monte-Carlo study for one parameter scenario and returns a
#' one-row data frame of operating characteristics, so rows for many scenarios
#' can simply be \code{rbind()}ed (e.g. over an \code{expand.grid()}). The
#' binomial analogue of \code{\link{SIM_norm}}. For each of the \code{n_sim}
#' replicates it samples \code{H} historical control studies and one current
#' trial from \code{\link{h1_binom}} (a beta-binomial random-effects model),
#' fits the borrowing analysis \code{\link{SCI_binom}}, and records the
#' per-contrast reject decisions, the prediction-interval decision and the
#' borrowing diagnostics; these are then averaged over replicates.
#'
#' @section Between-study variation:
#' The beta-binomial between-study variation is set with exactly one of
#' \code{rho}, \code{sd_h} or \code{var_h}. \code{rho} is the intra-class
#' correlation (the canonical parameter, invariant across \code{pi}). The
#' between-study variance of the true cluster proportions is
#' \eqn{\pi(1-\pi)\rho}, so a between-study standard deviation \code{sd_h} or
#' variance \code{var_h} is converted to \code{rho} via
#' \eqn{\rho = \mathrm{sd\_h}^2 / (\pi(1-\pi))} and
#' \eqn{\rho = \mathrm{var\_h} / (\pi(1-\pi))}, using the historical \code{pi}
#' (so the model holds \code{rho}, not the SD, constant across the historical
#' and current populations). This requires \eqn{\mathrm{var\_h} < \pi(1-\pi)}
#' (equivalently \code{rho < 1}).
#'
#' @param n_sim number of Monte-Carlo replicates (simulated trials).
#' @param H number of historical control studies.
#' @param n number of trials (patients) per arm (control and each treatment).
#' @param pi grand success probability of the data-generating process (shared by
#'   historical and current controls -- i.e. no prior-data conflict).
#' @param rho intra-class correlation of the beta-binomial model (between-study
#'   variation). Supply exactly one of \code{rho}, \code{sd_h}, \code{var_h}.
#' @param sd_h between-study standard deviation of the true cluster proportions
#'   (converted to \code{rho}; see the between-study variation section).
#' @param var_h between-study variance of the true cluster proportions
#'   (converted to \code{rho}).
#' @param shift true treatment-vs-control effect(s) on the scale given by
#'   \code{shift_scale}; \code{length(shift)} fixes the number of treatment
#'   groups. All-zero -> null (FWER); non-zero -> alternative (power).
#' @param shift_scale scale on which \code{shift} acts: \code{"prob"} (default,
#'   risk difference) or \code{"logit"} (log odds ratio). See \code{\link{h1_binom}}.
#' @param pi_cur grand success probability of the CURRENT trial. \code{NULL}
#'   (default) = \code{pi}, i.e. no prior-data conflict. A value != \code{pi}
#'   shifts the whole current trial (concurrent control AND every treatment arm)
#'   relative to the historical controls, inducing a prior-data conflict while
#'   leaving the treatment-vs-control effects unchanged.
#' @param prior_est analysis for each replicate: \code{"MAP"} (default) or
#'   \code{"EB"} borrowing (passed to \code{\link{SCI_binom}}), or
#'   \code{"ignore_HCD"} -- a non-borrowing frequentist competitor fitted here
#'   from the current trial alone (logistic GLM + unequal-variance Dunnett-type
#'   multiple comparison on odds ratios; requires \code{contr = "OR"}). Under
#'   \code{"ignore_HCD"} the borrowing-specific output columns (\code{pi_t1e},
#'   \code{pi_power}, \code{shrink_*}, \code{eff_n_*}, \code{post_w_*}, ...) are
#'   \code{NA}, so the runs \code{rbind()} together.
#' @param type multiple-comparison type passed through to \code{\link{SCI_binom}}:
#'   either \code{"Dunnett"} (many-to-one, each treatment vs. control; the
#'   default) or \code{"Tukey"} (all pairwise comparisons). Under \code{"Tukey"}
#'   the per-contrast \code{reject_*} output columns are keyed by the full
#'   pairwise contrast labels.
#' @param contr \code{"RR"} (ratio of proportions; default), \code{"OR"}
#'   (odds ratio) or \code{"pi_diff"} (difference of proportions).
#' @param limits \code{"two.sided"} (default), \code{"lower"} or
#'   \code{"upper"} (drives the reject decision; passed to \code{\link{SCI_binom}}).
#' @param base index of the control group in the contrast matrix.
#' @param robust weight of the weakly-informative (robust) prior component.
#' @param tau_prior,beta_prior optional overrides for the \code{gMAP} between-study
#'   heterogeneity and intercept priors, passed to \code{\link{SCI_binom}} (MAP
#'   path only). \code{NULL} (the default) uses the built-in heuristics.
#' @param alpha 1 - simultaneous confidence level.
#' @param n_mcmc number of Monte-Carlo draws inside each borrowing fit (passed to
#'   \code{\link{SCI_binom}}).
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
#'       least one shift != 0; \code{NA} under the global null.}
#'     \item{\code{pi_t1e}}{type-1 error of the prediction interval (rate the
#'       concurrent control proportion falls outside it), reported only under NO
#'       prior-data conflict (\code{pi_cur == pi}); \code{NA} otherwise.}
#'     \item{\code{pi_power}}{the same rate reported only WITH conflict
#'       (\code{pi_cur != pi}); \code{NA} under no conflict.}
#'     \item{\code{pi_reject_lower}, \code{pi_reject_upper}}{PI rejection rate
#'       split by border.}
#'     \item{\code{shrink_abs}, \code{shrink_rel}}{mean absolute / relative
#'       shrinkage of the concurrent control proportion towards the borrowed prior.}
#'     \item{\code{width_mean}, \code{eff_n_mean}, \code{eff_n_post_mean},
#'       \code{post_w_mean}, \code{post_w_median}}{further borrowing diagnostics
#'       averaged over replicates.}
#'     \item{\code{reject_*}, \code{reject_lower_*}, \code{reject_upper_*}}{
#'       per-contrast rejection rates.}
#'   }
#'
#' @seealso \code{\link{SCI_binom}}, \code{\link{h1_binom}}, \code{\link{SIM_norm}}
#'
#' @examples
#' \donttest{
#' # Null scenario -> FWER and PI type-1 error
#' SIM_binom(n_sim = 100, H = 5, n = 40, pi = 0.3, rho = 0.05,
#'           shift = c(0, 0), seed = 1, parallel = TRUE, n_cores = 2)
#'
#' # Prior-data conflict -> PI power to detect the shift
#' SIM_binom(n_sim = 100, H = 5, n = 40, pi = 0.3, rho = 0.05,
#'           shift = c(0, 0.15), pi_cur = 0.55, seed = 1,
#'           parallel = TRUE, n_cores = 2)
#'
#' # Non-borrowing frequentist competitor (logistic GLM + unequal-variance
#' # Dunnett on odds ratios; requires contr = "OR"). rbind() with a borrowing run
#' # to compare; the shared seed keeps the current trials identical (paired).
#' SIM_binom(n_sim = 100, H = 5, n = 40, pi = 0.3, rho = 0.05,
#'           shift = c(0, 0), prior_est = "ignore_HCD", contr = "OR",
#'           seed = 1, parallel = TRUE, n_cores = 2)
#' }
#' @export
SIM_binom <- function(n_sim,
                      H,
                      n,
                      pi,
                      rho         = NULL,
                      sd_h        = NULL,
                      var_h       = NULL,
                      shift       = c(0, 0),
                      shift_scale = c("prob", "logit"),
                      pi_cur      = NULL,
                      prior_est   = "MAP",
                      type        = "Dunnett",
                      contr       = "RR",
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

        if(!is.numeric(pi) || length(pi) != 1 || pi <= 0 || pi >= 1)
                stop("'pi' must be a single number in (0, 1).")
        if(!prior_est %in% c("MAP", "EB", "ignore_HCD"))
                stop("'prior_est' must be 'MAP', 'EB' or 'ignore_HCD'.")
        if(prior_est == "ignore_HCD" && contr != "OR")
                stop("prior_est = 'ignore_HCD' supports only contr = 'OR' for the binomial endpoint.")

        #-----------------------------------------------------------------------
        # Between-study variation: convert sd_h / var_h to the canonical rho.
        # var_between(pi_h) = pi (1 - pi) rho, so rho = var / (pi (1 - pi)); the
        # historical pi is used, keeping rho (not the SD) constant across the
        # historical and current populations.

        n_supplied <- sum(!is.null(rho), !is.null(sd_h), !is.null(var_h))
        if(n_supplied != 1)
                stop("Supply exactly one of 'rho', 'sd_h' or 'var_h'.")

        if(!is.null(sd_h))  rho <- sd_h^2 / (pi * (1 - pi))
        if(!is.null(var_h)) rho <- var_h  / (pi * (1 - pi))

        if(!is.numeric(rho) || length(rho) != 1 || rho <= 0 || rho >= 1){
                stop("The between-study variation is out of range: 'rho' must be in (0, 1), ",
                     "i.e. a between-study variance below pi*(1-pi) = ",
                     signif(pi * (1 - pi), 3), " (sd below ", signif(sqrt(pi * (1 - pi)), 3), ").")
        }

        # Between-study SD / variance of the true cluster proportions (reported)
        var_h_between <- pi * (1 - pi) * rho
        sd_h_between  <- sqrt(var_h_between)

        # Current-trial grand probability. NULL -> pi (no prior-data conflict).
        if(is.null(pi_cur)) pi_cur <- pi
        if(!is.numeric(pi_cur) || length(pi_cur) != 1 || pi_cur <= 0 || pi_cur >= 1)
                stop("'pi_cur' must be a single number in (0, 1).")

        # Force all arguments before the parallel dispatch, so run_one()'s closure
        # captures their VALUES rather than promises bound to the caller's
        # environment. Without this, PSOCK workers cannot resolve an argument
        # passed as an unevaluated variable expression (e.g. H = H).
        for(.nm in names(formals())) force(get(.nm, envir = environment(), inherits = FALSE))

        # Silence RBesT / rstan chatter; genuine errors still propagate.
        quiet <- function(expr) suppressWarnings(suppressMessages(expr))

        #-----------------------------------------------------------------------
        # Per-replicate operating characteristics from one fitted SCI_binom() object

        op_char <- function(fit){

                ccg_mean  <- fit$newdat$mean[1]
                pw        <- fit$posterior["w", ]
                pa        <- fit$posterior["a", ]
                pb        <- fit$posterior["b", ]
                post_mean <- sum(pw * pa / (pa + pb))

                # Per-contrast reject decisions, resolved by credible limit
                reject_lower <- fit$comp$lower.ci > null_val
                reject_upper <- fit$comp$upper.ci < null_val
                width        <- fit$comp$upper.ci - fit$comp$lower.ci

                # Prediction-interval decision for the concurrent control proportion
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
        # One Monte-Carlo replicate: sample HCD + current trial, fit SCI_binom(),
        # reduce to operating characteristics. Returns NULL if the fit fails.

        run_one <- function(stream){

                if(!is.null(stream)){
                        RNGkind("L'Ecuyer-CMRG")
                        assign(".Random.seed", stream, envir = .GlobalEnv)
                }

                # Historical control studies (control arm only, beta-binomial)
                hist_raw <- h1_binom(pi = pi, H = H, n = n, rho = rho,
                                     shift = NULL, shift_scale = shift_scale)
                histdat  <- data.frame(events     = hist_raw$events,
                                       non_events = hist_raw$n - hist_raw$events)

                # Current trial (one study). pi_cur shifts it vs. the HCD (conflict);
                # the treatment arms share the drawn cluster proportion.
                cur_raw <- h1_binom(pi = pi_cur, H = 1, n = n, rho = rho,
                                    shift = shift, shift_scale = shift_scale)
                newdat  <- data.frame(group      = as.character(cur_raw$arm),
                                      events     = cur_raw$events,
                                      non_events = cur_raw$n - cur_raw$events,
                                      stringsAsFactors = FALSE)

                # Non-borrowing frequentist competitor: logistic GLM + unequal-
                # variance Dunnett-type multiple comparison (odds ratios), current
                # trial only. Borrowing-specific quantities are NA.
                if(prior_est == "ignore_HCD"){
                        fc <- try(quiet(.freq_competitor(newdat, distr = "betabinomial", type = type,
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

                fit <- try(quiet(SCI_binom(histdat = histdat, newdat = newdat,
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

                # analysis packages on every worker (SCI_binom() / h1_binom() need these)
                parallel::clusterEvalQ(cl, suppressWarnings(suppressPackageStartupMessages({
                        library(RBesT); library(multcomp)
                        library(mratios); library(BSagri); library(predint)
                })))
                # mirror the caller's RBesT MCMC options; run Stan chains
                # sequentially inside each worker to avoid oversubscription
                rbest_opts <- options()[grep("^RBesT\\.", names(options()))]
                parallel::clusterCall(cl, function(o){ options(o); options(mc.cores = 1L) },
                                      rbest_opts)
                # ship the analysis functions (all scalars/params travel with
                # run_one()'s closure environment)
                parallel::clusterExport(cl, varlist = c("SCI_binom", "h1_binom", ".freq_competitor"),
                                        envir = environment())

                oc_list <- parallel::parLapply(cl, streams, run_one)
        } else {
                oc_list <- lapply(streams, run_one)
        }

        # Drop replicates whose borrowing fit failed (SCI_binom() wrapped in try())
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

        # Prediction-interval type-1 error (and per border). Under NO conflict
        # (pi_cur == pi) this is the PI's type-1 error; WITH conflict it is its
        # power to detect the whole-trial shift -- mirroring the fwer / power split.
        pi_reject       <- mean(vapply(oc_list, function(r) r$reject_pi, logical(1)))
        pi_reject_lower <- mean(vapply(oc_list, function(r) r$reject_pl, logical(1)))
        pi_reject_upper <- mean(vapply(oc_list, function(r) r$reject_pu, logical(1)))

        pi_conflict <- !isTRUE(all.equal(pi_cur, pi))          # whole-trial shift vs HCD
        pi_t1e   <- if(!pi_conflict) pi_reject else NA_real_   # no conflict -> type-1 error
        pi_power <- if( pi_conflict) pi_reject else NA_real_   # conflict    -> power to detect it

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
                          pi              = pi,
                          pi_cur          = pi_cur,
                          rho             = rho,
                          sd_h_between    = sd_h_between,
                          var_h_between   = var_h_between,
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

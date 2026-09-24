#' Plot a MultiMAP object
#'
#' Single S3 plot method for MultiMAP objects, for normal
#' (\code{\link{SCI_norm}} / \code{\link{PI_norm}}), binomial
#' (\code{\link{SCI_binom}} / \code{\link{PI_binom}}) and Poisson
#' (\code{\link{SCI_pois}}) endpoints; the graphic is chosen via \code{which},
#' depending on how the object was produced:
#' \itemize{
#'   \item \code{PI} object (class \code{"PI"}): only the prediction-interval
#'     plot is defined, so \code{which} is ignored (a value other than 1 warns).
#'   \item \code{SCI} object (class \code{"SCI"}):
#'     \code{which = 1} prediction-interval plot (default);
#'     \code{which = 2} simultaneous credible-interval plot;
#'     \code{which = 3} the prior, likelihood and posterior of the concurrent
#'     control overlaid (visualises shrinkage).
#' }
#' For a normal endpoint the prediction-interval plot is a continuous density
#' with the interval shaded; for a binomial endpoint it is the discrete
#' beta-binomial predictive drawn on the proportion scale as a bar chart, and for
#' a Poisson endpoint the discrete gamma-Poisson predictive drawn on the rate
#' scale -- in both discrete cases the bars inside the interval are highlighted.
#'
#' @param x a MultiMAP object from \code{\link{SCI_norm}}, \code{\link{PI_norm}},
#'   \code{\link{SCI_binom}} or \code{\link{PI_binom}}.
#' @param which \code{1} (default) prediction-interval plot, \code{2}
#'   simultaneous-CI plot, \code{3} prior/likelihood/posterior plot (both
#'   \code{2} and \code{3} are for \code{SCI_norm} objects only).
#' @param arm for multi-arm objects (\code{histdat} supplied as a named list)
#'   only: which borrowed arm the per-arm prediction-interval (\code{which = 1})
#'   and prior/likelihood/posterior (\code{which = 3}) plots refer to. \code{NULL}
#'   (the default) draws every borrowed arm as its own sub-panel
#'   (\code{ggplot2::facet_wrap}); giving an arm name restricts the plot to that
#'   single arm. Ignored for single-arm objects and for the simultaneous-CI plot
#'   (\code{which = 2}).
#' @param ... forwarded to the prediction-interval plot, namely:
#'   \describe{
#'     \item{\code{shade_pi}}{\code{TRUE} (default) shades the prediction interval
#'       \code{x$pred_int} (lower ... upper) as a blue area under the density
#'       curve; \code{FALSE} omits it. One-sided intervals (an \code{NA} border)
#'       are shaded up to the corresponding edge of the plotting range.}
#'     \item{\code{show_ccg}}{\code{TRUE} (default) draws the concurrent control
#'       mean (\code{x$newdat$mean[1]}) as a vertical line -- black if covered by
#'       the prediction interval, red if not -- for \code{SCI_norm} input;
#'       \code{FALSE} omits it. Silently skipped for \code{PI_norm} input.}
#'     \item{\code{prob}}{central probability mass of the mixture used for the
#'       x-range (as in \code{RBesT::plot.mix}; default 0.99).}
#'     \item{\code{fill}, \code{PI_alpha}}{fill colour and opacity of the shaded
#'       prediction interval.}
#'   }
#'
#' @return A ggplot object.
#'
#' @examples
#' \donttest{
#' histdat <- data.frame(mean = c(9.8, 10.2, 9.5, 10.0),
#'                       sd   = c(2.1, 1.9, 2.3, 2.0),
#'                       n    = c(40, 55, 38, 47))
#' newdat  <- data.frame(group = c("control", "d1", "d454"),
#'                       mean  = c(10.1, 11.4, 12.7),
#'                       sd    = c(2.0, 2.2, 2.1),
#'                       n     = c(30, 30, 30))
#' res <- SCI_norm(histdat, newdat, contr = "mean_diff", seed = 1)
#' plot(res, which = 1)   # prediction interval
#' plot(res, which = 2)   # simultaneous credible intervals
#' plot(res, which = 3)   # prior / likelihood / posterior
#' }
#' @export
plot.MultiMAP <- function(x, which = 1, arm = NULL, ...){

        #-------------------- helper: prediction-interval plot -----------------
        # Normal endpoint: density of the (continuous) prediction distribution
        # x$pred_dist as a black line, with the prediction interval x$pred_int
        # shaded as a blue area under the curve.
        # Binomial endpoint: the prediction distribution is the discrete
        # beta-binomial predictive for the number of events; it is drawn on the
        # proportion scale (events / n) as a bar chart, with the bars inside the
        # prediction interval highlighted.
        # For SCI_*() input the observed concurrent control mean / proportion is
        # added as a vertical line, black if covered by the prediction interval
        # and red if not.
        plot_PI <- function(x,
                            shade_pi = TRUE,
                            show_ccg = TRUE,
                            prob     = 0.99,
                            fill     = "lightblue",
                            PI_alpha = 0.55){

                mix <- x$pred_dist
                pi  <- x$pred_int

                #============ binomial: discrete beta-binomial predictive =======
                if(identical(x$distr, "betabinomial")){

                        # Number of future patients: the largest count in the
                        # predictive support (its 100% quantile).
                        n0 <- round(unname(RBesT::qmix(mix, 1)))
                        rr    <- 0:n0
                        df    <- data.frame(prop = rr / n0,
                                            pmf  = RBesT::dmix(mix, rr))

                        # Which bars fall inside the prediction interval
                        # (an NA border counts as unbounded).
                        lo <- if(is.na(pi["lower"])) -Inf else unname(pi["lower"])
                        hi <- if(is.na(pi["upper"]))  Inf else unname(pi["upper"])
                        df$inside <- df$prop >= lo & df$prop <= hi

                        # Focus the x axis on the informative region (the mass, the
                        # prediction interval and the observed value); drop the
                        # empty tail a heavy-tailed predictive would otherwise show.
                        obs <- if(inherits(x, "SCI")) x$newdat$mean[1] else NA
                        win <- .info_window(df$prop, df$pmf, include = c(lo, hi, obs))
                        df  <- df[df$prop >= win[1] & df$prop <= win[2], , drop = FALSE]

                        g <- ggplot(df, aes(x = prop, y = pmf)) + theme_bw()

                        if(isTRUE(shade_pi)){
                                g <- g + geom_col(aes(fill = inside), width = 0.9 / n0) +
                                        scale_fill_manual(values = c(`TRUE` = fill, `FALSE` = "grey80"),
                                                          guide = "none")
                        } else {
                                g <- g + geom_col(width = 0.9 / n0, fill = "grey60")
                        }

                        if(isTRUE(show_ccg) && inherits(x, "SCI")){
                                ccg_mean <- x$newdat$mean[1]
                                covered  <- (is.na(pi["lower"]) || ccg_mean >= pi["lower"]) &&
                                            (is.na(pi["upper"]) || ccg_mean <= pi["upper"])
                                g <- g + geom_vline(xintercept = ccg_mean,
                                                    colour     = if(isTRUE(unname(covered))) "black" else "red")
                        }

                        return(g + ylab("probability") + xlab("proportion"))
                }

                #============ poisson: discrete gamma-Poisson predictive ========
                if(identical(x$distr, "poisson")){

                        # Control exposure the predictive was built for (used to
                        # move from event counts to the rate scale).
                        exposure <- if(inherits(x, "SCI")) x$newdat[1, 3] else 1

                        # Count support up to a high predictive quantile
                        rr <- 0:ceiling(unname(RBesT::qmix(mix, 0.9999)))
                        df <- data.frame(rate = rr / exposure,
                                         pmf  = RBesT::dmix(mix, rr))

                        lo <- if(is.na(pi["lower"])) -Inf else unname(pi["lower"])
                        hi <- if(is.na(pi["upper"]))  Inf else unname(pi["upper"])
                        df$inside <- df$rate >= lo & df$rate <= hi

                        # Focus the x axis on the informative region (see above).
                        obs <- if(inherits(x, "SCI")) x$newdat$mean[1] else NA
                        win <- .info_window(df$rate, df$pmf, include = c(lo, hi, obs))
                        df  <- df[df$rate >= win[1] & df$rate <= win[2], , drop = FALSE]

                        g <- ggplot(df, aes(x = rate, y = pmf)) + theme_bw()

                        if(isTRUE(shade_pi)){
                                g <- g + geom_col(aes(fill = inside), width = 0.9 / exposure) +
                                        scale_fill_manual(values = c(`TRUE` = fill, `FALSE` = "grey80"),
                                                          guide = "none")
                        } else {
                                g <- g + geom_col(width = 0.9 / exposure, fill = "grey60")
                        }

                        if(isTRUE(show_ccg) && inherits(x, "SCI")){
                                ccg_mean <- x$newdat$mean[1]
                                covered  <- (is.na(pi["lower"]) || ccg_mean >= pi["lower"]) &&
                                            (is.na(pi["upper"]) || ccg_mean <= pi["upper"])
                                g <- g + geom_vline(xintercept = ccg_mean,
                                                    colour     = if(isTRUE(unname(covered))) "black" else "red")
                        }

                        return(g + ylab("probability") + xlab("rate"))
                }

                #============ normal: continuous predictive density =============

                # x-range: central 'prob' mass of the mixture (as RBesT::plot.mix)
                plow     <- (1 - prob) / 2
                interval <- RBesT::qmix(mix, c(plow, 1 - plow))

                # density of the full predictive mixture (the black curve)
                dens_mix <- function(z) RBesT::dmix(mix, z)

                g <- ggplot(data.frame(x = interval), aes(x = x)) + theme_bw()

                #---- shaded prediction interval (behind the density line) -----
                if(isTRUE(shade_pi)){
                        lo <- if(is.na(pi["lower"])) interval[1] else unname(pi["lower"])
                        hi <- if(is.na(pi["upper"])) interval[2] else unname(pi["upper"])
                        xx <- seq(lo, hi, length.out = 501)
                        g <- g + geom_area(data = data.frame(x = xx, y = RBesT::dmix(mix, xx)),
                                           mapping = aes(x = x, y = y),
                                           fill = fill, alpha = PI_alpha,
                                           inherit.aes = FALSE)
                }

                #---- black mixture density ------------------------------------
                g <- g + stat_function(geom = "line", fun = dens_mix,
                                       n = 501)

                #---- concurrent control mean (SCI_norm() input only) ----------
                # Black if covered by the prediction interval, red otherwise
                # (a border of NA counts as unbounded).
                if(isTRUE(show_ccg) && inherits(x, "SCI")){
                        ccg_mean <- x$newdat$mean[1]
                        covered  <- (is.na(pi["lower"]) || ccg_mean >= pi["lower"]) &&
                                    (is.na(pi["upper"]) || ccg_mean <= pi["upper"])
                        g <- g + geom_vline(xintercept = ccg_mean,
                                            colour     = if(isTRUE(unname(covered))) "black" else "red")
                }

                g + ylab("density") + xlab("parameter")
        }

        #-------------------- helper: simultaneous-CI plot ---------------------
        # Point estimates + simultaneous credible intervals of the contrasts,
        # with the null value marked (0 for difference-type contrasts, 1 for
        # ratio-type contrasts). Works for both normal endpoints (mean_diff /
        # mean_ratio) and binomial endpoints (pi_diff / RR / OR); the null
        # value and title are derived from x$contr via the shared helpers.
        plot_SCI <- function(x){

                null_val <- .null_value(x$contr)   # 0 for differences, 1 for ratios
                ttl      <- .contr_title(x$contr)

                sci <- get_output(x, which = "SCI")

                # Order the y-axis to match the row order of the results ($comp):
                # ggplot draws the first factor level at the bottom, so reverse the
                # levels to put the first comparison at the top (as in the table).
                sci$comp <- factor(sci$comp, levels = rev(as.character(sci$comp)))

                # One-sided intervals carry an infinite border; drop it so the
                # point range is drawn open-ended rather than clipping the scale.
                sci$lower.ci[!is.finite(sci$lower.ci)] <- NA
                sci$upper.ci[!is.finite(sci$upper.ci)] <- NA

                ggplot(sci, aes(y = comp)) +
                        theme_bw() +
                        geom_pointrange(aes(x = estimate, xmin = lower.ci, xmax = upper.ci)) +
                        geom_vline(xintercept = null_val, linetype = "dashed") +
                        ylab("Comparison") + xlab("SCI") + ggtitle(ttl)
        }

        #-------------------- helper: prior, likelihood and posterior ----------
        # Overlay, on the concurrent-control-mean scale, the three densities:
        #   * Likelihood  -- the concurrent control data alone, N(ybar, se),
        #   * Prior       -- the robustified MAP prior actually updated (prior_rob),
        #   * Posterior   -- the borrowed posterior of the control mean.
        # The posterior sitting between likelihood and prior visualises shrinkage;
        # dashed vertical lines mark the raw (likelihood) and posterior means so
        # the shift is easy to read off.
        plot_prior_lik_post <- function(x){

                prior_mix <- x$prior_rob    # operative (robustified) prior
                post_mix  <- x$posterior

                lik_lab   <- "Likelihood (CCG data)"
                prior_lab <- "Prior (robustified)"
                post_lab  <- "Posterior"

                cols <- c("#E69F00", "#56B4E9", "#009E73")
                names(cols) <- c(lik_lab, prior_lab, post_lab)

                #============ binomial: densities on the proportion scale =======
                if(identical(x$distr, "betabinomial")){

                        # Concurrent control data: events r out of n.
                        r0    <- x$newdat[1, 2]
                        n0    <- x$newdat[1, 2] + x$newdat[1, 3]
                        ccg_mean <- x$newdat$mean[1]

                        # exact posterior-mixture mean (beta comp. mean = a/(a+b))
                        post_mean <- sum(post_mix["w", ] *
                                                 post_mix["a", ] / (post_mix["a", ] + post_mix["b", ]))

                        # Normalised likelihood of the proportion, i.e. the flat
                        # (Beta(1, 1)) prior updated with the CCG data:
                        #   L(p) proportional to Beta(r + 1, n - r + 1).
                        # Evaluate the three densities on a fine [0, 1] grid, then
                        # focus the x axis on the informative region (where their
                        # envelope is >= 1% of its peak, widened to the observed and
                        # posterior means). This drops the empty tail a heavy-tailed
                        # MAP prior would otherwise stretch across [0, 1].
                        xx0   <- seq(0, 1, length.out = 1024)
                        d_lik <- stats::dbeta(xx0, r0 + 1, n0 - r0 + 1)
                        d_pri <- RBesT::dmix(prior_mix, xx0)
                        d_pos <- RBesT::dmix(post_mix,  xx0)
                        win   <- .info_window(xx0, pmax(d_lik, d_pri, d_pos),
                                              include = c(ccg_mean, post_mean))
                        keep  <- xx0 >= win[1] & xx0 <= win[2]
                        dens  <- rbind(
                                data.frame(x = xx0[keep], density = d_lik[keep], Distribution = lik_lab),
                                data.frame(x = xx0[keep], density = d_pri[keep], Distribution = prior_lab),
                                data.frame(x = xx0[keep], density = d_pos[keep], Distribution = post_lab))
                        dens$Distribution <- factor(dens$Distribution,
                                levels = c(lik_lab, prior_lab, post_lab))

                        return(ggplot(dens, aes(x = x, y = density, colour = Distribution)) +
                                       theme_bw() +
                                       geom_line(linewidth = 1) +
                                       geom_vline(xintercept = ccg_mean,  colour = cols[[lik_lab]],
                                                  linetype = "dashed") +
                                       geom_vline(xintercept = post_mean, colour = cols[[post_lab]],
                                                  linetype = "dashed") +
                                       scale_colour_manual(values = cols) +
                                       ylab("density") + xlab("concurrent control proportion") +
                                       ggtitle("Prior, likelihood and posterior"))
                }

                #============ poisson: densities on the rate scale ==============
                if(identical(x$distr, "poisson")){

                        # Concurrent control data: y_c events over exposure t_c.
                        y_ccg    <- x$newdat[1, 2]
                        t_ccg    <- x$newdat[1, 3]
                        ccg_mean <- x$newdat$mean[1]

                        # exact posterior-mixture mean (gamma comp. mean = a/b)
                        post_mean <- sum(post_mix["w", ] * post_mix["a", ] / post_mix["b", ])

                        # Normalised likelihood of the rate, i.e. the flat prior
                        # updated with the CCG data: L(lambda) proportional to
                        # Gamma(y_c + 1, rate = t_c). Evaluate the three densities on
                        # a fine grid over a generous range, then focus the x axis on
                        # the informative region (envelope >= 1% of its peak, widened
                        # to the observed and posterior means).
                        hi0 <- max(RBesT::qmix(prior_mix, 0.999),
                                   RBesT::qmix(post_mix,  0.999),
                                   stats::qgamma(0.999, shape = y_ccg + 1, rate = t_ccg))
                        xx0   <- seq(0, hi0, length.out = 1024)
                        d_lik <- dgamma(xx0, shape = y_ccg + 1, rate = t_ccg)
                        d_pri <- RBesT::dmix(prior_mix, xx0)
                        d_pos <- RBesT::dmix(post_mix,  xx0)
                        win   <- .info_window(xx0, pmax(d_lik, d_pri, d_pos),
                                              include = c(ccg_mean, post_mean))
                        keep  <- xx0 >= max(0, win[1]) & xx0 <= win[2]
                        dens  <- rbind(
                                data.frame(x = xx0[keep], density = d_lik[keep], Distribution = lik_lab),
                                data.frame(x = xx0[keep], density = d_pri[keep], Distribution = prior_lab),
                                data.frame(x = xx0[keep], density = d_pos[keep], Distribution = post_lab))
                        dens$Distribution <- factor(dens$Distribution,
                                levels = c(lik_lab, prior_lab, post_lab))

                        return(ggplot(dens, aes(x = x, y = density, colour = Distribution)) +
                                       theme_bw() +
                                       geom_line(linewidth = 1) +
                                       geom_vline(xintercept = ccg_mean,  colour = cols[[lik_lab]],
                                                  linetype = "dashed") +
                                       geom_vline(xintercept = post_mean, colour = cols[[post_lab]],
                                                  linetype = "dashed") +
                                       scale_colour_manual(values = cols) +
                                       ylab("density") + xlab("concurrent control rate") +
                                       ggtitle("Prior, likelihood and posterior"))
                }

                #============ normal: densities on the mean scale ===============

                # concurrent control likelihood: N(ybar, se)
                ccg_mean <- x$newdat$mean[1]
                ccg_se   <- x$newdat$sd[1] / sqrt(x$newdat$n[1])

                # x-range covering all three densities
                rng <- range(RBesT::qmix(prior_mix, c(0.005, 0.995)),
                             RBesT::qmix(post_mix,  c(0.005, 0.995)),
                             ccg_mean + c(-4, 4) * ccg_se)
                xx  <- seq(rng[1], rng[2], length.out = 512)

                dens <- rbind(
                        data.frame(x = xx, density = dnorm(xx, ccg_mean, ccg_se),
                                   Distribution = lik_lab),
                        data.frame(x = xx, density = RBesT::dmix(prior_mix, xx),
                                   Distribution = prior_lab),
                        data.frame(x = xx, density = RBesT::dmix(post_mix, xx),
                                   Distribution = post_lab))
                dens$Distribution <- factor(dens$Distribution,
                        levels = c(lik_lab, prior_lab, post_lab))

                # exact posterior-mixture mean (raw mean = ccg_mean); gap = shrinkage
                post_mean <- sum(post_mix["w", ] * post_mix["m", ])

                ggplot(dens, aes(x = x, y = density, colour = Distribution)) +
                        theme_bw() +
                        geom_line(linewidth = 1) +
                        geom_vline(xintercept = ccg_mean,  colour = cols[[lik_lab]],
                                   linetype = "dashed") +
                        geom_vline(xintercept = post_mean, colour = cols[[post_lab]],
                                   linetype = "dashed") +
                        scale_colour_manual(values = cols) +
                        ylab("density") + xlab("concurrent control mean") +
                        ggtitle("Prior, likelihood and posterior")
        }

        #-------------------- helper: multi-arm prediction-interval facets -----
        # Same graphic as plot_PI() but with one sub-panel per borrowed arm
        # (facet_wrap, free scales). Discrete endpoints (binomial / Poisson) are
        # bar charts on the proportion / rate scale; the normal endpoint is a
        # continuous density with the prediction interval shaded. For SCI input
        # the observed value is added per panel (black if covered, red if not).
        plot_PI_facets <- function(x,
                                   shade_pi = TRUE,
                                   show_ccg = TRUE,
                                   prob     = 0.99,
                                   fill     = "lightblue",
                                   PI_alpha = 0.55){

                arms     <- x$borrow_arms
                from_SCI <- inherits(x, "SCI")
                grp      <- if(from_SCI) as.character(x$newdat$group) else NULL
                bars     <- x$distr %in% c("betabinomial", "poisson")

                curve_l <- list(); shade_l <- list(); vline_l <- list()

                for(a in arms){
                        mix <- x$pred_dist[[a]]
                        pi  <- x$pred_int[[a]]
                        idx <- if(from_SCI) match(a, grp) else NA
                        obs <- if(from_SCI) x$newdat$mean[idx] else NA
                        lo  <- if(is.na(pi["lower"])) -Inf else unname(pi["lower"])
                        hi  <- if(is.na(pi["upper"]))  Inf else unname(pi["upper"])

                        if(identical(x$distr, "betabinomial")){
                                n0 <- round(unname(RBesT::qmix(mix, 1)))
                                rr <- 0:n0
                                df <- data.frame(x = rr / n0, y = RBesT::dmix(mix, rr),
                                                 width = 0.9 / n0)
                        } else if(identical(x$distr, "poisson")){
                                exposure <- if(from_SCI) x$newdat[idx, 3] else 1
                                rr <- 0:ceiling(unname(RBesT::qmix(mix, 0.9999)))
                                df <- data.frame(x = rr / exposure, y = RBesT::dmix(mix, rr),
                                                 width = 0.9 / exposure)
                        } else {
                                plow     <- (1 - prob) / 2
                                interval <- RBesT::qmix(mix, c(plow, 1 - plow))
                                xx <- seq(interval[1], interval[2], length.out = 501)
                                df <- data.frame(x = xx, y = RBesT::dmix(mix, xx), width = NA_real_)
                        }
                        df$inside <- df$x >= lo & df$x <= hi

                        # Discrete endpoints: focus each panel's x axis on the
                        # informative region (the mass, the prediction interval and
                        # the observed value); drop the empty heavy tail.
                        if(bars){
                                win <- .info_window(df$x, df$y, include = c(lo, hi, obs))
                                df  <- df[df$x >= win[1] & df$x <= win[2], , drop = FALSE]
                        }
                        df$arm    <- a
                        curve_l[[a]] <- df

                        # continuous endpoint: shaded prediction interval per panel
                        if(!bars && isTRUE(shade_pi)){
                                slo <- if(is.infinite(lo)) min(df$x) else lo
                                shi <- if(is.infinite(hi)) max(df$x) else hi
                                xx2 <- seq(slo, shi, length.out = 501)
                                shade_l[[a]] <- data.frame(x = xx2, y = RBesT::dmix(mix, xx2), arm = a)
                        }
                        # observed concurrent value per panel (SCI input only)
                        if(isTRUE(show_ccg) && from_SCI){
                                covered <- (is.infinite(lo) || obs >= lo) &&
                                           (is.infinite(hi) || obs <= hi)
                                vline_l[[a]] <- data.frame(arm = a, xintercept = obs,
                                                           col = if(isTRUE(unname(covered))) "black" else "red")
                        }
                }

                curve      <- do.call(rbind, curve_l)
                curve$arm  <- factor(curve$arm, levels = arms)
                xlab_txt   <- switch(x$distr, betabinomial = "proportion",
                                     poisson = "rate", "parameter")

                g <- ggplot() + theme_bw()

                if(bars){
                        if(isTRUE(shade_pi)){
                                g <- g + geom_col(data = curve,
                                                  aes(x = x, y = y, fill = inside, width = width)) +
                                        scale_fill_manual(values = c(`TRUE` = fill, `FALSE` = "grey80"),
                                                          guide = "none")
                        } else {
                                g <- g + geom_col(data = curve, aes(x = x, y = y, width = width),
                                                  fill = "grey60")
                        }
                        ylab_txt <- "probability"
                } else {
                        if(length(shade_l)){
                                shade     <- do.call(rbind, shade_l)
                                shade$arm <- factor(shade$arm, levels = arms)
                                g <- g + geom_area(data = shade, aes(x = x, y = y),
                                                   fill = fill, alpha = PI_alpha)
                        }
                        g <- g + geom_line(data = curve, aes(x = x, y = y))
                        ylab_txt <- "density"
                }

                if(length(vline_l)){
                        vdf     <- do.call(rbind, vline_l)
                        vdf$arm <- factor(vdf$arm, levels = arms)
                        g <- g + geom_vline(data = vdf, aes(xintercept = xintercept, colour = col)) +
                                scale_colour_identity(guide = "none")
                }

                g + facet_wrap(~ arm, scales = "free") + ylab(ylab_txt) + xlab(xlab_txt)
        }

        #-------------------- helper: multi-arm prior/likelihood/posterior facets
        # Same overlay as plot_prior_lik_post() but with one sub-panel per
        # borrowed arm (SCI input only; needs the per-arm posterior + newdat).
        plot_plp_facets <- function(x){

                arms      <- x$borrow_arms
                lik_lab   <- "Likelihood (CCG data)"
                prior_lab <- "Prior (robustified)"
                post_lab  <- "Posterior"
                cols        <- c("#E69F00", "#56B4E9", "#009E73")
                names(cols) <- c(lik_lab, prior_lab, post_lab)
                grp <- as.character(x$newdat$group)

                dens_l <- list(); ccg_l <- list(); post_l <- list()

                for(a in arms){
                        prior_mix <- x$prior_rob[[a]]
                        post_mix  <- x$posterior[[a]]
                        idx       <- match(a, grp)

                        if(identical(x$distr, "betabinomial")){
                                r0 <- x$newdat[idx, 2]; n0 <- x$newdat[idx, 2] + x$newdat[idx, 3]
                                ccg_mean <- x$newdat$mean[idx]
                                post_mean <- sum(post_mix["w", ] *
                                                 post_mix["a", ] / (post_mix["a", ] + post_mix["b", ]))
                                # Fine [0, 1] grid, then focus on the informative
                                # region (envelope >= 1% of its peak, widened to the
                                # observed and posterior means).
                                xx0   <- seq(0, 1, length.out = 1024)
                                d_lik <- stats::dbeta(xx0, r0 + 1, n0 - r0 + 1)
                                d_pri <- RBesT::dmix(prior_mix, xx0)
                                d_pos <- RBesT::dmix(post_mix,  xx0)
                                win   <- .info_window(xx0, pmax(d_lik, d_pri, d_pos),
                                                      include = c(ccg_mean, post_mean))
                                keep  <- xx0 >= win[1] & xx0 <= win[2]
                                dens  <- rbind(
                                        data.frame(x = xx0[keep], density = d_lik[keep], Distribution = lik_lab),
                                        data.frame(x = xx0[keep], density = d_pri[keep], Distribution = prior_lab),
                                        data.frame(x = xx0[keep], density = d_pos[keep], Distribution = post_lab))
                        } else if(identical(x$distr, "poisson")){
                                y_ccg <- x$newdat[idx, 2]; t_ccg <- x$newdat[idx, 3]
                                ccg_mean <- x$newdat$mean[idx]
                                post_mean <- sum(post_mix["w", ] * post_mix["a", ] / post_mix["b", ])
                                hi0 <- max(RBesT::qmix(prior_mix, 0.999),
                                           RBesT::qmix(post_mix,  0.999),
                                           stats::qgamma(0.999, shape = y_ccg + 1, rate = t_ccg))
                                xx0   <- seq(0, hi0, length.out = 1024)
                                d_lik <- dgamma(xx0, shape = y_ccg + 1, rate = t_ccg)
                                d_pri <- RBesT::dmix(prior_mix, xx0)
                                d_pos <- RBesT::dmix(post_mix,  xx0)
                                win   <- .info_window(xx0, pmax(d_lik, d_pri, d_pos),
                                                      include = c(ccg_mean, post_mean))
                                keep  <- xx0 >= max(0, win[1]) & xx0 <= win[2]
                                dens  <- rbind(
                                        data.frame(x = xx0[keep], density = d_lik[keep], Distribution = lik_lab),
                                        data.frame(x = xx0[keep], density = d_pri[keep], Distribution = prior_lab),
                                        data.frame(x = xx0[keep], density = d_pos[keep], Distribution = post_lab))
                        } else {
                                ccg_mean <- x$newdat$mean[idx]
                                ccg_se   <- x$newdat$sd[idx] / sqrt(x$newdat$n[idx])
                                rng <- range(RBesT::qmix(prior_mix, c(0.005, 0.995)),
                                             RBesT::qmix(post_mix,  c(0.005, 0.995)),
                                             ccg_mean + c(-4, 4) * ccg_se)
                                xx  <- seq(rng[1], rng[2], length.out = 512)
                                dens <- rbind(
                                        data.frame(x = xx, density = dnorm(xx, ccg_mean, ccg_se),
                                                   Distribution = lik_lab),
                                        data.frame(x = xx, density = RBesT::dmix(prior_mix, xx),
                                                   Distribution = prior_lab),
                                        data.frame(x = xx, density = RBesT::dmix(post_mix, xx),
                                                   Distribution = post_lab))
                                post_mean <- sum(post_mix["w", ] * post_mix["m", ])
                        }

                        dens$arm     <- a
                        dens_l[[a]]  <- dens
                        ccg_l[[a]]   <- data.frame(arm = a, xintercept = ccg_mean)
                        post_l[[a]]  <- data.frame(arm = a, xintercept = post_mean)
                }

                dens              <- do.call(rbind, dens_l)
                dens$Distribution <- factor(dens$Distribution, levels = c(lik_lab, prior_lab, post_lab))
                dens$arm          <- factor(dens$arm, levels = arms)
                ccg_df            <- do.call(rbind, ccg_l);  ccg_df$arm  <- factor(ccg_df$arm,  levels = arms)
                post_df           <- do.call(rbind, post_l); post_df$arm <- factor(post_df$arm, levels = arms)

                xlab_txt <- switch(x$distr, betabinomial = "concurrent control proportion",
                                   poisson = "concurrent control rate", "concurrent control mean")

                ggplot(dens, aes(x = x, y = density, colour = Distribution)) +
                        theme_bw() +
                        geom_line(linewidth = 1) +
                        geom_vline(data = ccg_df,  aes(xintercept = xintercept),
                                   colour = cols[[lik_lab]],  linetype = "dashed") +
                        geom_vline(data = post_df, aes(xintercept = xintercept),
                                   colour = cols[[post_lab]], linetype = "dashed") +
                        scale_colour_manual(values = cols) +
                        facet_wrap(~ arm, scales = "free") +
                        ylab("density") + xlab(xlab_txt) +
                        ggtitle("Prior, likelihood and posterior")
        }

        #-------------------- dispatch -----------------------------------------
        # SCI_norm() objects carry class "SCI"; PI_norm() objects carry class "PI".
        from_SCI <- inherits(x, "SCI")

        # Multi-arm objects hold the per-arm prior/PI/posterior components as named
        # lists. The SCI plot (which = 2) is a single table and is unaffected; the
        # prediction-interval (1) and prior/likelihood/posterior (3) plots are
        # per-arm. By default (arm = NULL) they show one sub-panel per borrowed arm
        # (facet_wrap); passing 'arm' instead selects a single arm and reshapes 'x'
        # to a single-arm view of it.
        if(!is.null(x$borrow_arms)){
                which_i <- if(from_SCI) as.integer(which)[1] else 1L
                per_arm <- which_i %in% c(1L, 3L)

                if(per_arm && is.null(arm)){
                        # one sub-panel per borrowed arm
                        if(which_i == 1L) return(plot_PI_facets(x, ...))
                        if(which_i == 3L) return(plot_plp_facets(x))
                }

                if(per_arm){
                        # a single, named arm: reshape 'x' to a single-arm view
                        a <- as.character(arm)[1]
                        if(!a %in% x$borrow_arms)
                                stop("'arm' must be one of the borrowed arms: ",
                                     paste(x$borrow_arms, collapse = ", "), ".")
                        x$pred_dist <- x$pred_dist[[a]]
                        x$pred_int  <- x$pred_int[[a]]
                        if(is.list(x$prior))     x$prior     <- x$prior[[a]]
                        if(is.list(x$prior_rob)) x$prior_rob <- x$prior_rob[[a]]
                        if(is.list(x$posterior)) x$posterior <- x$posterior[[a]]
                        # Put the selected arm's current data in row 1 (used for the
                        # likelihood curve and the observed-value line). PI_*()
                        # objects carry no newdat, so only SCI objects need this.
                        if(from_SCI)
                                x$newdat <- x$newdat[match(a, as.character(x$newdat$group)), , drop = FALSE]
                }
                # which = 2 falls through unchanged (arm-independent SCI table plot)
        }

        if(!from_SCI){
                # PI_norm() object -> only the prediction-interval plot exists
                if(!identical(as.integer(which), 1L))
                        warning("PI_norm() objects only support the prediction-interval plot (which = 1); ignoring 'which'.")
                return(plot_PI(x, ...))
        }

        # SCI_norm() object -> choose the graphic via 'which'
        which <- as.integer(which)[1]
        if(which == 1L) return(plot_PI(x, ...))     # prediction-interval plot
        if(which == 2L) return(plot_SCI(x))         # simultaneous CI plot
        if(which == 3L) return(plot_prior_lik_post(x))   # prior, likelihood and posterior

        stop("'which' must be 1 (prediction interval), 2 (simultaneous CI) or 3 (prior, likelihood and posterior).")
}

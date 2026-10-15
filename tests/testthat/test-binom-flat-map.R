# Regression (issue H1): the default logit-scale HalfNormal tau heuristic
# round(1 / sqrt(p_bar (1 - p_bar)), 1) / 2 grew without bound as p_bar -> 0.
# For all-zero / rare historical data it produced a tau of ~6-11, an absurdly
# diffuse (near-flat, bimodal) MAP prior: prior mean ~0.40 from 0/150 events,
# a prediction interval spanning [0, 1], ~20x-inflated odds ratios in SCI_binom,
# and a crash inside RBesT::qmix for rare non-zero data. The default is now a
# fixed HalfNormal(0, 1) on the logit scale, with a guardrail warning when the
# predictive still spans almost all of [0, 1]. These use the MAP path (gMAP),
# so they skip on CRAN and use a small n_mcmc.

test_that("SCI_pois-style blow-up is gone: all-zero HCD gives a sane MAP prior", {
        skip_on_cran()
        hd <- data.frame(events = c(0, 0, 0), non_events = c(50, 50, 50))
        p  <- suppressMessages(suppressWarnings(
                PI_binom(hd, n0 = 50, prior_est = "MAP", seed = 1)))

        # default tau scale is now the fixed 1 (was 6.1 for 3 zero studies)
        expect_equal(as.numeric(p$tau_prior$prior), 1)
        # prior is concentrated near 0, not pulled to ~0.40
        expect_lt(as.numeric(summary(p$prior)["mean"]), 0.1)
        # the predictive interval must not span almost all of [0, 1]
        span <- as.numeric(p$pred_int["upper"]) - as.numeric(p$pred_int["lower"])
        expect_lt(span, 0.9)
})

test_that("rare non-zero HCD no longer crashes PI_binom", {
        skip_on_cran()
        hd2 <- data.frame(events = c(1, 0, 0, 1, 0), non_events = c(49, 50, 50, 49, 50))
        expect_error(
                suppressMessages(suppressWarnings(
                        PI_binom(hd2, n0 = 50, prior_est = "MAP", seed = 1))),
                NA)
})

test_that("SCI_binom odds ratio is not inflated for a zero-count control", {
        skip_on_cran()
        hd <- data.frame(events = c(0, 0, 0), non_events = c(50, 40, 45))
        nd <- data.frame(group = c("C", "T2"), events = c(0, 3), non_events = c(50, 47))
        cc <- suppressMessages(suppressWarnings(
                SCI_binom(hd, nd, contr = "OR", n_mcmc = 2000, seed = 1)))$comp
        # was ~326 with upper CI ~1.45e5 under the old flat-prior heuristic; with
        # the H1 tau fix and the M7 data-centred intercept it is ~87 [4.7, 4320]
        # (a genuinely large OR for a 0-event control vs a 3-event arm, but far
        # from the degenerate value). Thresholds stay tight enough to catch a
        # regression to the flat-prior behaviour.
        expect_lt(cc$estimate[1], 150)
        expect_lt(cc$upper.ci[1], 1e4)
})

test_that("guardrail warns when a forced wide tau makes the predictive near-flat", {
        skip_on_cran()
        # Capture + filter warnings so unrelated gMAP/automixfit warnings neither
        # fail the check nor leak out of the test.
        warns_of <- function(expr) {
                w <- character(0)
                withCallingHandlers(
                        suppressMessages(expr),
                        warning = function(cnd) { w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning") })
                w
        }

        hd <- data.frame(events = c(0, 0, 0), non_events = c(50, 50, 50))
        w_flat <- warns_of(PI_binom(hd, n0 = 50, prior_est = "MAP", seed = 1, tau_prior = 6))
        expect_true(any(grepl("spans almost all", w_flat)))

        # a well-identified data set must NOT trigger this warning
        hw <- data.frame(events = c(10, 12, 9, 11, 8), non_events = c(40, 38, 41, 39, 42))
        w_ok <- warns_of(PI_binom(hw, n0 = 50, prior_est = "MAP", seed = 1))
        expect_false(any(grepl("spans almost all", w_ok)))
})

# Regression (issue C1): non-borrowed Poisson arms must not collapse a zero-count
# arm to rate 0. The old near-flat Gamma(0.001, 0.001) prior put ~49% of a
# zero-count arm's posterior draws at exactly 0, giving spurious significant rate
# ratios (FWER ~0.42 under the null). The weakly-informative Gamma(0.5, eps) prior
# fixes it. These use the EB path (no gMAP), so they run fast.

test_that("SCI_pois: a zero-count non-borrowed arm does not collapse to rate 0", {
        hp <- data.frame(events = c(2, 1, 3, 2, 1, 2, 3, 2), offset = 5)
        np <- data.frame(group  = c("0", "10", "50"), events = c(2, 0, 4), offset = 5)
        r  <- suppressWarnings(SCI_pois(hp, np, prior_est = "EB", seed = 1, n_mcmc = 5000))

        # the zero-count arm (column 2 of simpost) must not pile mass at exactly 0
        expect_lt(mean(r$simpost[, 2] == 0), 0.01)   # was ~0.49 before the fix

        # its rate-ratio SCI (10 / 0) must cover the null value 1 (not spurious)
        rr <- r$comp[r$comp$comp == "10/0", ]
        expect_true(is.finite(rr$lower.ci) && is.finite(rr$upper.ci))
        expect_lte(rr$lower.ci, 1)
        expect_gte(rr$upper.ci, 1)
})

test_that("SCI_pois: all-zero HCD + a zero-count control gives finite bounds", {
        hz <- data.frame(events = rep(0, 8), offset = 5)
        nz <- data.frame(group = c("0", "10"), events = c(0, 2), offset = 5)
        r  <- suppressWarnings(SCI_pois(hz, nz, prior_est = "EB", seed = 1, n_mcmc = 5000))
        expect_true(all(is.finite(c(r$comp$lower.ci, r$comp$upper.ci))))
})

test_that("SCI_pois multi-arm: a non-borrowed zero-count arm stays finite", {
        hl <- list("0" = data.frame(events = c(2, 1, 3, 2), offset = 5))
        nm <- data.frame(group = c("0", "10", "50"), events = c(2, 0, 4), offset = 5)
        r  <- suppressWarnings(SCI_pois(hl, nm, prior_est = "EB", seed = 1, n_mcmc = 5000))
        rr <- r$comp[r$comp$comp == "10/0", ]
        expect_true(is.finite(rr$lower.ci) && is.finite(rr$upper.ci))
        expect_lte(rr$lower.ci, 1)
        expect_gte(rr$upper.ci, 1)
})

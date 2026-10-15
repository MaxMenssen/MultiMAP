# Regression (issue M4): a non-integer 'base' passed validation (only
# 1 <= base <= nrow(newdat) was checked) and silently produced contrasts with no
# control coefficient (for difference/OR via multcomp::contrMat) or a truncated
# base (for RR/rate_ratio via mratios::contrMatRatio). 'base' must now be a whole
# number. EB path -> fast, no gMAP.

test_that("SCI_* reject a non-integer base and accept an integer one", {
        hb <- data.frame(events = c(2, 2, 1, 5, 3, 5), non_events = c(48, 47, 49, 44, 51, 44))
        nb <- data.frame(group = c("0", "10", "50"), events = c(3, 4, 8), non_events = c(44, 36, 39))
        expect_error(SCI_binom(hb, nb, prior_est = "EB", base = 1.5, n_mcmc = 1000, seed = 1),
                     "whole-number")
        expect_error(suppressWarnings(SCI_binom(hb, nb, prior_est = "EB", base = 2, n_mcmc = 1000, seed = 1)), NA)

        hp <- data.frame(events = c(10, 12, 8, 15, 11), offset = 5)
        np <- data.frame(group = c("0", "100", "200"), events = c(10, 15, 20), offset = 5)
        expect_error(SCI_pois(hp, np, prior_est = "EB", contr = "rate_diff",
                              base = 1.5, n_mcmc = 1000, seed = 1), "whole-number")

        hn <- data.frame(mean = c(10, 11, 9), sd = c(2, 2, 2), n = c(30, 30, 30))
        nn <- data.frame(group = c("0", "1"), mean = c(10, 12), sd = c(2, 2), n = c(30, 30))
        expect_error(SCI_norm(hn, nn, base = 1.5, n_mcmc = 1000, seed = 1), "whole-number")
})

test_that("SCI_binom multi-arm also rejects a non-integer base", {
        hb <- data.frame(events = c(2, 2, 1, 5, 3, 5), non_events = c(48, 47, 49, 44, 51, 44))
        nb <- data.frame(group = c("0", "10"), events = c(3, 4), non_events = c(44, 36))
        expect_error(SCI_binom(list("0" = hb), nb, prior_est = "EB",
                               base = 1.5, n_mcmc = 1000, seed = 1), "whole-number")
})

test_that("SIM_* reject a non-integer base early (before running replicates)", {
        expect_error(SIM_binom(n_sim = 4, H = 5, n = 40, pi = 0.3, rho = 0.05,
                               shift = c(0, 0), prior_est = "EB", base = 1.5,
                               parallel = FALSE, seed = 1), "whole-number")
        expect_error(SIM_pois(n_sim = 4, H = 6, offset = 3, lambda = 8, shape = 12,
                              shift = c(0, 0), prior_est = "EB", base = 1.5,
                              parallel = FALSE, seed = 1), "whole-number")
})

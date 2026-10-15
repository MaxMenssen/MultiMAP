# Regression (issue M5): tibbles pass the is.data.frame() validation but their
# `[`, `[,` indexing returns a 1-column tibble rather than a vector, so the
# single-arm SCI_binom / SCI_pois / PI_binom paths crashed with cryptic base /
# RBesT errors (the multi-arm paths happened to work). Inputs are now coerced to
# a base data frame at entry, so tibbles give the same results as data frames.
# EB path -> fast, no gMAP.

test_that("tibble inputs give the same results as data.frame inputs (single-arm)", {
        skip_if_not_installed("tibble")
        tb <- tibble::as_tibble

        hb <- data.frame(events = c(2, 2, 1, 5, 3, 5, 2, 2),
                         non_events = c(48, 47, 49, 44, 51, 44, 54, 53))
        nb <- data.frame(group = c("0", "10", "50"),
                         events = c(3, 4, 8), non_events = c(44, 36, 39))

        # SCI_binom: tibble newdat and tibble histdat both crashed before
        s_df <- suppressWarnings(SCI_binom(hb, nb, prior_est = "EB", n_mcmc = 1500, seed = 1))
        s_tb <- suppressWarnings(SCI_binom(tb(hb), tb(nb), prior_est = "EB", n_mcmc = 1500, seed = 1))
        expect_equal(s_tb$comp, s_df$comp, tolerance = 1e-9)
        expect_s3_class(s_tb$newdat, "data.frame")        # coerced to a base data frame

        # PI_binom: tibble histdat crashed before
        p_df <- suppressWarnings(PI_binom(hb, n0 = 50, prior_est = "EB", n_mcmc = 1500, seed = 1))
        p_tb <- suppressWarnings(PI_binom(tb(hb), n0 = 50, prior_est = "EB", n_mcmc = 1500, seed = 1))
        expect_equal(p_tb$pred_int, p_df$pred_int, tolerance = 1e-9)

        # SCI_pois: tibble newdat crashed before
        hp <- data.frame(events = c(5, 7, 4, 6, 8, 5, 6, 7), offset = c(10, 12, 9, 11, 13, 10, 11, 12))
        np <- data.frame(group = c("0", "10", "50"), events = c(6, 8, 9), offset = rep(10, 3))
        sp_df <- suppressWarnings(SCI_pois(hp, np, prior_est = "EB", n_mcmc = 1500, seed = 1))
        sp_tb <- suppressWarnings(SCI_pois(hp, tb(np), prior_est = "EB", n_mcmc = 1500, seed = 1))
        expect_equal(sp_tb$comp, sp_df$comp, tolerance = 1e-9)
})

test_that("tibble inputs also work in the multi-arm path", {
        skip_if_not_installed("tibble")
        tb <- tibble::as_tibble
        hb <- data.frame(events = c(2, 2, 1, 5, 3, 5, 2, 2),
                         non_events = c(48, 47, 49, 44, 51, 44, 54, 53))
        nb <- data.frame(group = c("0", "10"), events = c(3, 4), non_events = c(44, 36))
        expect_error(
                suppressWarnings(SCI_binom(list("0" = tb(hb)), tb(nb),
                                           prior_est = "EB", n_mcmc = 1500, seed = 1)), NA)
})

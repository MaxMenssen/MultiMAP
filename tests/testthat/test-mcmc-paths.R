# MAP / normal paths run gMAP() MCMC and are therefore slow: skip on CRAN and
# whenever RBesT is unavailable. They still run in local / CI checks.

test_that("SCI_norm (MAP) produces a normal SCI object", {
        skip_on_cran()
        skip_if_not_installed("RBesT")

        histdat <- data.frame(mean = c(9.8, 10.2, 9.5, 10.0),
                              sd   = c(2.1, 1.9, 2.3, 2.0),
                              n    = c(40, 55, 38, 47))
        newdat  <- data.frame(group = c("control", "d1", "d454"),
                              mean  = c(10.1, 11.4, 12.7),
                              sd    = c(2.0, 2.2, 2.1),
                              n     = c(30, 30, 30))
        res <- suppressWarnings(SCI_norm(histdat, newdat, contr = "mean_diff", seed = 1))

        expect_true(inherits(res, "SCI") && inherits(res, "MultiMAP_norm"))
        expect_identical(res$distr, "normal")
        expect_equal(nrow(res$comp), 2)
        expect_s3_class(plot(res, which = 3), "ggplot")
})

test_that("SCI_norm validates the contrast type", {
        skip_on_cran()
        skip_if_not_installed("RBesT")
        histdat <- data.frame(mean = c(9.8, 10.2, 9.5), sd = c(2.1, 1.9, 2.3),
                              n = c(40, 55, 38))
        newdat  <- data.frame(group = c("control", "d1"), mean = c(10.1, 11.4),
                              sd = c(2.0, 2.2), n = c(30, 30))
        # match.arg() rejects an unknown contrast before any MCMC runs
        expect_error(SCI_norm(histdat, newdat, contr = "mean_prop"))
})

test_that("SCI_binom (MAP) produces a binomial SCI object", {
        skip_on_cran()
        skip_if_not_installed("RBesT")

        histdat <- data.frame(events     = c(11, 8, 14, 9),
                              non_events = c(29, 32, 26, 31))
        newdat  <- data.frame(group      = c("control", "d1", "d454"),
                              events     = c(8, 14, 17),
                              non_events = c(22, 16, 13))
        res <- suppressWarnings(SCI_binom(histdat, newdat, prior_est = "MAP",
                                          contr = "pi_diff", seed = 1))
        expect_true(inherits(res, "SCI") && inherits(res, "MultiMAP_binom"))
        expect_equal(nrow(res$comp), 2)
})

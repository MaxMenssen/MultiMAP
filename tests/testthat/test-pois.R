# The empirical-Bayes (EB) Poisson path needs no gMAP MCMC, so these run fast.
# Overdispersed historical counts keep the moment estimator away from its floor.

make_pois_data <- function(){
        list(histdat = data.frame(events = c(4, 22, 6, 35, 9, 28, 3, 26),
                                  offset = rep(2, 8)),
             newdat  = data.frame(group  = c("control", "d1", "d454"),
                                  events = c(10, 18, 21),
                                  offset = c(3, 3, 3)))
}

test_that("SCI_pois (EB) returns a well-formed MultiMAP object", {
        d   <- make_pois_data()
        res <- suppressWarnings(SCI_pois(d$histdat, d$newdat, prior_est = "EB",
                                         contr = "rate_diff", seed = 1))

        expect_s3_class(res, "MultiMAP")
        expect_true(inherits(res, "SCI"))
        expect_true(inherits(res, "MultiMAP_pois"))
        expect_identical(res$distr, "poisson")

        expect_named(res$comp, c("comp", "estimate", "lower.ci", "upper.ci"))
        expect_equal(nrow(res$comp), 2)
        expect_named(res$pred_int, c("median", "lower", "upper"))
})

test_that("get_output('ess') returns the three Poisson effective sample sizes", {
        d   <- make_pois_data()
        res <- suppressWarnings(SCI_pois(d$histdat, d$newdat, prior_est = "EB", seed = 1))
        ess <- get_output(res, which = "ess")
        expect_named(ess, c("MAP ESS", "Robust. MAP ESS", "Post. ESS"))
        expect_length(ess, 3)
})

test_that("SCI_pois validates the contrast type", {
        d <- make_pois_data()
        expect_error(SCI_pois(d$histdat, d$newdat, prior_est = "EB", contr = "prop"))
})

test_that("SCI_pois rejects an unknown prior_est", {
        d <- make_pois_data()
        expect_error(SCI_pois(d$histdat, d$newdat, prior_est = "ML"))
})

test_that("shared S3 methods work on a Poisson (EB) fit", {
        d   <- make_pois_data()
        res <- suppressWarnings(SCI_pois(d$histdat, d$newdat, prior_est = "EB",
                                         contr = "rate_ratio", seed = 1))
        expect_s3_class(plot(res, which = 1), "ggplot")
        expect_s3_class(plot(res, which = 2), "ggplot")
        expect_s3_class(plot(res, which = 3), "ggplot")
        expect_output(summary(res), "rate")
        expect_output(summary(res), "ratio of rates")
})

test_that("SCI_pois (MAP) produces a Poisson SCI object", {
        skip_on_cran()
        skip_if_not_installed("RBesT")
        d   <- make_pois_data()
        res <- suppressWarnings(SCI_pois(d$histdat, d$newdat, prior_est = "MAP",
                                         contr = "rate_diff", seed = 1))
        expect_true(inherits(res, "SCI") && inherits(res, "MultiMAP_pois"))
        expect_equal(nrow(res$comp), 2)
})

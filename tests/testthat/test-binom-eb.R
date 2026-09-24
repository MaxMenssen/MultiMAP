# The empirical-Bayes (EB) prior needs no gMAP MCMC, so these tests run fast
# and are safe to run everywhere (including CRAN).

make_binom_data <- function(){
        list(histdat = data.frame(events     = c(11, 8, 14, 9),
                                  non_events = c(29, 32, 26, 31)),
             newdat  = data.frame(group      = c("control", "d1", "d454"),
                                  events     = c(8, 14, 17),
                                  non_events = c(22, 16, 13)))
}

test_that("SCI_binom (EB) returns a well-formed MultiMAP object", {
        d   <- make_binom_data()
        res <- suppressWarnings(SCI_binom(d$histdat, d$newdat, prior_est = "EB",
                                          contr = "pi_diff", seed = 1))

        expect_s3_class(res, "MultiMAP")
        expect_true(inherits(res, "SCI"))
        expect_true(inherits(res, "MultiMAP_binom"))
        expect_identical(res$distr, "betabinomial")

        # comp table: one row per treatment, standardised columns
        expect_named(res$comp, c("comp", "estimate", "lower.ci", "upper.ci"))
        expect_equal(nrow(res$comp), 2)

        # prediction interval on the proportion scale
        expect_named(res$pred_int, c("median", "lower", "upper"))
        expect_true(res$pred_int["lower"] >= 0 && res$pred_int["upper"] <= 1)
})

test_that("get_output('ess') returns the three effective sample sizes", {
        d   <- make_binom_data()
        res <- suppressWarnings(SCI_binom(d$histdat, d$newdat, prior_est = "EB", seed = 1))
        ess <- get_output(res, which = "ess")
        expect_named(ess, c("MAP ESS", "Robust. MAP ESS", "Post. ESS"))
        expect_length(ess, 3)
        expect_true(all(ess >= 0))
})

test_that("SCI_binom validates the contrast type", {
        d <- make_binom_data()
        expect_error(SCI_binom(d$histdat, d$newdat, prior_est = "EB",
                               contr = "pi_ratio"))   # not a valid contr
})

test_that("SCI_binom rejects an unknown prior_est (e.g. the old 'Lui2000')", {
        d <- make_binom_data()
        expect_error(SCI_binom(d$histdat, d$newdat, prior_est = "Lui2000"))
})

test_that("PI_binom (EB) returns a PI object with an NA posterior", {
        d      <- make_binom_data()
        m_pi   <- suppressWarnings(PI_binom(d$histdat, n0 = 30, prior_est = "EB", seed = 1))
        expect_true(inherits(m_pi, "PI"))
        expect_true(inherits(m_pi, "MultiMAP_binom"))
        expect_named(m_pi$pred_int, c("median", "lower", "upper"))
        expect_true(is.na(m_pi$eff_n_post))
})

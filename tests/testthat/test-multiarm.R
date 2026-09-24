# Multi-arm borrowing: 'histdat' supplied as a named list borrows for more than
# one arm. The binomial / Poisson EB priors need no gMAP MCMC, so those tests run
# fast and everywhere; the normal endpoint (MAP only) is guarded with
# skip_on_cran().

# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------

make_binom_ma <- function(){
        list(hist = list(control = data.frame(events     = c(11, 8, 14, 9),
                                              non_events = c(29, 32, 26, 31)),
                         d454    = data.frame(events     = c(18, 21, 16, 20),
                                              non_events = c(22, 19, 24, 20))),
             newdat = data.frame(group      = c("control", "d1", "d454"),
                                 events     = c(8, 14, 17),
                                 non_events = c(22, 16, 13)))
}

make_pois_ma <- function(){
        list(hist = list(control = data.frame(events = c(18, 22, 15, 25, 20), offset = 5),
                         d2      = data.frame(events = c(10, 13, 9, 12, 11),  offset = 5)),
             newdat = data.frame(group  = c("control", "d1", "d2"),
                                 events = c(95, 60, 62),
                                 offset = 5))
}

# ---------------------------------------------------------------------------
# Single-arm behaviour is preserved when histdat is a plain data frame
# ---------------------------------------------------------------------------

test_that("a data-frame histdat keeps the single-arm object (no borrow_arms)", {
        d   <- make_binom_ma()
        res <- suppressWarnings(SCI_binom(d$hist$control, d$newdat, prior_est = "EB", seed = 1))
        expect_null(res$borrow_arms)
        expect_length(res$eff_n, 1L)
        expect_named(res$pred_int, c("median", "lower", "upper"))
})

# ---------------------------------------------------------------------------
# Binomial multi-arm (EB)
# ---------------------------------------------------------------------------

test_that("SCI_binom multi-arm borrows the named arms and keeps per-arm slots", {
        d   <- make_binom_ma()
        res <- suppressWarnings(SCI_binom(d$hist, d$newdat, prior_est = "EB",
                                          contr = "RR", seed = 1))

        expect_s3_class(res, "MultiMAP")
        expect_true(inherits(res, "SCI"))
        expect_identical(res$distr, "betabinomial")

        # both named arms are borrowed, in newdat order
        expect_identical(res$borrow_arms, c("control", "d454"))

        # per-arm effective sample sizes are a named vector keyed by arm
        expect_named(res$eff_n, c("control", "d454"))
        expect_length(res$eff_n, 2L)

        # per-arm prediction intervals are a named list keyed by arm
        expect_type(res$pred_int, "list")
        expect_named(res$pred_int, c("control", "d454"))
        expect_named(res$pred_int$control, c("median", "lower", "upper"))

        # priors / posteriors are per-arm named lists
        expect_named(res$prior,     c("control", "d454"))
        expect_named(res$posterior, c("control", "d454"))

        # contrasts: one row per treatment, standardised columns
        expect_named(res$comp, c("comp", "estimate", "lower.ci", "upper.ci"))
        expect_equal(nrow(res$comp), 2)
})

test_that("get_output('ess') is a per-arm data frame in multi-arm mode", {
        d   <- make_binom_ma()
        res <- suppressWarnings(SCI_binom(d$hist, d$newdat, prior_est = "EB", seed = 1))
        ess <- get_output(res, which = "ess")
        expect_s3_class(ess, "data.frame")
        expect_identical(rownames(ess), c("control", "d454"))
        expect_named(ess, c("MAP ESS", "Robust. MAP ESS", "Post. ESS"))
})

test_that("summary() and print() work on a multi-arm object", {
        d   <- make_binom_ma()
        res <- suppressWarnings(SCI_binom(d$hist, d$newdat, prior_est = "EB", seed = 1))
        expect_output(summary(res), "Multi-arm borrowing for 2 arm")
        expect_output(print(res),   "control")
})

test_that("multi-arm plot draws one sub-panel per borrowed arm", {
        skip_if_not_installed("ggplot2")
        d   <- make_binom_ma()
        res <- suppressWarnings(SCI_binom(d$hist, d$newdat, prior_est = "EB", seed = 1))

        n_panels <- function(g) nrow(ggplot2::ggplot_build(g)$layout$layout)

        # which = 1 (prediction interval) and which = 3 (prior/lik/post): one
        # panel per borrowed arm by default (the bug fix).
        expect_equal(n_panels(plot(res, which = 1)), 2L)
        expect_equal(n_panels(plot(res, which = 3)), 2L)

        # a named 'arm' restricts to that single arm ...
        expect_equal(n_panels(plot(res, which = 1, arm = "d454")), 1L)
        # ... and the simultaneous-CI plot (which = 2) is a single panel.
        expect_equal(n_panels(plot(res, which = 2)), 1L)
})

test_that("PI_binom multi-arm returns one prediction interval per arm", {
        d    <- make_binom_ma()
        m_pi <- suppressWarnings(PI_binom(d$hist, n0 = 40, prior_est = "EB", seed = 1))
        expect_true(inherits(m_pi, "PI"))
        expect_identical(m_pi$borrow_arms, c("control", "d454"))
        expect_named(m_pi$pred_int, c("control", "d454"))
        expect_true(all(is.na(m_pi$eff_n_post)))
})

test_that("PI_binom multi-arm accepts a named n0 vector", {
        d    <- make_binom_ma()
        m_pi <- suppressWarnings(PI_binom(d$hist, n0 = c(control = 40, d454 = 80),
                                          prior_est = "EB", seed = 1))
        expect_named(m_pi$pred_int, c("control", "d454"))
})

# ---------------------------------------------------------------------------
# Poisson multi-arm (EB)
# ---------------------------------------------------------------------------

test_that("SCI_pois multi-arm borrows the named arms and keeps per-arm slots", {
        d   <- make_pois_ma()
        res <- suppressWarnings(SCI_pois(d$hist, d$newdat, prior_est = "EB",
                                         contr = "rate_ratio", seed = 1))

        expect_true(inherits(res, "SCI"))
        expect_identical(res$distr, "poisson")
        expect_identical(res$borrow_arms, c("control", "d2"))
        expect_named(res$eff_n, c("control", "d2"))
        expect_type(res$pred_int, "list")
        expect_named(res$pred_int, c("control", "d2"))
        expect_equal(nrow(res$comp), 2)
})

test_that("PI_pois multi-arm returns one prediction interval per arm", {
        d    <- make_pois_ma()
        m_pi <- suppressWarnings(PI_pois(d$hist, offset0 = 5, prior_est = "EB", seed = 1))
        expect_true(inherits(m_pi, "PI"))
        expect_identical(m_pi$borrow_arms, c("control", "d2"))
        expect_named(m_pi$pred_int, c("control", "d2"))
})

# ---------------------------------------------------------------------------
# Validation of the named-list histdat
# ---------------------------------------------------------------------------

test_that("multi-arm histdat is validated against newdat's groups", {
        d <- make_binom_ma()

        # an unnamed list is rejected
        unnamed <- unname(d$hist)
        expect_error(SCI_binom(unnamed, d$newdat, prior_est = "EB"),
                     "named")

        # a name that matches no arm is rejected
        bad <- d$hist
        names(bad) <- c("control", "nope")
        expect_error(SCI_binom(bad, d$newdat, prior_est = "EB"),
                     "do not match")
})

test_that("a non-borrowed arm is allowed (uses the flat-prior posterior)", {
        d   <- make_binom_ma()
        # borrow only the control arm; d454 is a treatment arm without HCD
        res <- suppressWarnings(SCI_binom(d$hist["control"], d$newdat,
                                          prior_est = "EB", contr = "RR", seed = 1))
        expect_identical(res$borrow_arms, "control")
        expect_length(res$eff_n, 1L)
        expect_equal(nrow(res$comp), 2)
})

# ---------------------------------------------------------------------------
# Normal multi-arm (MAP only -> needs MCMC)
# ---------------------------------------------------------------------------

test_that("SCI_norm multi-arm produces a per-arm MultiMAP object", {
        skip_on_cran()
        skip_if_not_installed("RBesT")

        hist <- list(control = data.frame(mean = c(9.8, 10.2, 9.5, 10.0),
                                          sd   = c(2.1, 1.9, 2.3, 2.0),
                                          n    = c(40, 55, 38, 47)),
                     d2      = data.frame(mean = c(12.1, 11.8, 12.4),
                                          sd   = c(2.0, 2.1, 1.9),
                                          n    = c(35, 40, 30)))
        newdat <- data.frame(group = c("control", "d1", "d2"),
                             mean  = c(10.1, 11.4, 12.2),
                             sd    = c(2.0, 2.2, 2.1),
                             n     = c(30, 30, 30))

        res <- suppressWarnings(SCI_norm(hist, newdat, contr = "mean_diff", seed = 1))
        expect_true(inherits(res, "SCI"))
        expect_identical(res$distr, "normal")
        expect_identical(res$borrow_arms, c("control", "d2"))
        expect_named(res$eff_n, c("control", "d2"))
        expect_named(res$pred_int, c("control", "d2"))
        expect_equal(nrow(res$comp), 2)
})

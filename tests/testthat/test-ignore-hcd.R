# The non-borrowing frequentist competitor lives only in the SIM_*() functions
# (prior_est = "ignore_HCD"): standard GLM + multcomp Dunnett on the current
# trial. No MCMC, so these run fast.

test_that("SCI_*/PI_* do NOT accept 'ignore_HCD' (competitor is SIM-only)", {
        hn <- data.frame(mean = c(9.8, 10.2, 9.5), sd = c(2.1, 1.9, 2.3), n = c(40, 55, 38))
        nn <- data.frame(group = c("control", "d1"), mean = c(10.1, 11.4), sd = c(2, 2.2), n = c(30, 30))
        expect_error(SCI_norm(hn, nn, prior_est = "ignore_HCD"))

        hb <- data.frame(events = c(11, 8, 14, 9), non_events = c(29, 32, 26, 31))
        expect_error(PI_binom(hb, n0 = 30, prior_est = "ignore_HCD"))
})

test_that("SIM_norm ignore_HCD returns a competitor row with NA borrowing columns", {
        r <- suppressWarnings(SIM_norm(n_sim = 30, H = 5, n = 30, mu = 10, sd_h = 2,
                                       sd_hi_mean = 2, sd_hi_cv = 0.5, shift = c(0, 0),
                                       prior_est = "ignore_HCD", seed = 1, parallel = FALSE))
        expect_equal(r$prior_est, "ignore_HCD")
        expect_false(is.na(r$fwer))          # computable
        expect_false(is.na(r$width_mean))    # computable (frequentist CI width)
        expect_true(is.na(r$pi_t1e))         # borrowing-specific -> NA
        expect_true(is.na(r$shrink_abs))
        expect_true(is.na(r$eff_n_mean))
        expect_true(is.na(r$post_w_mean))
})

test_that("SIM_binom / SIM_pois ignore_HCD work on their natural contrast", {
        rb <- suppressWarnings(SIM_binom(n_sim = 30, H = 6, n = 40, pi = 0.3, rho = 0.05,
                                         shift = c(0, 0), prior_est = "ignore_HCD",
                                         contr = "OR", seed = 1, parallel = FALSE))
        expect_equal(rb$prior_est, "ignore_HCD")
        expect_false(is.na(rb$fwer))
        expect_true(is.na(rb$pi_power) || is.na(rb$pi_t1e))

        rp <- suppressWarnings(SIM_pois(n_sim = 30, H = 6, offset = 3, lambda = 8, shape = 12,
                                        shift = c(0, 0), prior_est = "ignore_HCD",
                                        contr = "rate_ratio", seed = 1, parallel = FALSE))
        expect_equal(rp$prior_est, "ignore_HCD")
        expect_false(is.na(rp$fwer))
})

test_that("ignore_HCD hard-stops on a non-natural contrast", {
        expect_error(SIM_binom(n_sim = 2, H = 6, n = 40, pi = 0.3, rho = 0.05,
                               shift = c(0, 0), prior_est = "ignore_HCD",
                               contr = "RR", parallel = FALSE))
        expect_error(SIM_pois(n_sim = 2, H = 6, offset = 3, lambda = 8, shape = 12,
                              shift = c(0, 0), prior_est = "ignore_HCD",
                              contr = "rate_diff", parallel = FALSE))
})

test_that("a borrowing run and an ignore_HCD run rbind() (same columns)", {
        r_eb  <- suppressWarnings(SIM_pois(n_sim = 6, H = 6, offset = 3, lambda = 8, shape = 12,
                                           shift = c(0, 0), prior_est = "EB", seed = 1, parallel = FALSE))
        r_ign <- suppressWarnings(SIM_pois(n_sim = 6, H = 6, offset = 3, lambda = 8, shape = 12,
                                           shift = c(0, 0), prior_est = "ignore_HCD", seed = 1, parallel = FALSE))
        expect_identical(names(r_eb), names(r_ign))
        expect_equal(nrow(rbind(r_eb, r_ign)), 2)
})

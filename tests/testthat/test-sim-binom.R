# SIM_binom with the EB prior needs no gMAP MCMC, so a small run is fast.

test_that("SIM_binom (EB) returns a one-row data frame with the key columns", {
        out <- SIM_binom(n_sim = 8, H = 5, n = 40, pi = 0.3, rho = 0.05,
                         shift = c(0, 0), prior_est = "EB", contr = "pi_diff",
                         seed = 1, parallel = FALSE)
        expect_s3_class(out, "data.frame")
        expect_equal(nrow(out), 1)
        expect_true(all(c("pi", "pi_cur", "rho", "sd_h_between", "var_h_between",
                          "fwer", "power", "pi_t1e", "pi_power") %in% names(out)))
        # Null scenario: FWER defined, power NA
        expect_false(is.na(out$fwer))
        expect_true(is.na(out$power))
})

test_that("SIM_binom converts sd_h / var_h to the same rho", {
        pi  <- 0.3
        rho <- 0.05
        sdi <- sqrt(pi * (1 - pi) * rho)
        o1 <- SIM_binom(n_sim = 3, H = 4, n = 40, pi = pi, sd_h  = sdi,
                        shift = c(0, 0), prior_est = "EB", seed = 1, parallel = FALSE)
        o2 <- SIM_binom(n_sim = 3, H = 4, n = 40, pi = pi, var_h = pi*(1-pi)*rho,
                        shift = c(0, 0), prior_est = "EB", seed = 1, parallel = FALSE)
        expect_equal(o1$rho, rho, tolerance = 1e-8)
        expect_equal(o2$rho, rho, tolerance = 1e-8)
})

test_that("SIM_binom requires exactly one between-study specification", {
        expect_error(SIM_binom(n_sim = 2, H = 3, n = 40, pi = 0.3,
                               shift = c(0, 0), prior_est = "EB", parallel = FALSE))
        expect_error(SIM_binom(n_sim = 2, H = 3, n = 40, pi = 0.3,
                               rho = 0.05, sd_h = 0.1, shift = c(0, 0),
                               prior_est = "EB", parallel = FALSE))
})

test_that("SIM_binom rejects a between-study variance that exceeds pi(1-pi)", {
        expect_error(SIM_binom(n_sim = 2, H = 3, n = 40, pi = 0.3,
                               var_h = 0.5, shift = c(0, 0), prior_est = "EB",
                               parallel = FALSE))
})

test_that("SIM_binom conflict scenario reports pi_power, not pi_t1e", {
        out <- SIM_binom(n_sim = 8, H = 5, n = 40, pi = 0.3, rho = 0.05,
                         shift = c(0, 0.15), pi_cur = 0.55, prior_est = "EB",
                         contr = "pi_diff", seed = 1, parallel = FALSE)
        expect_true(is.na(out$pi_t1e))
        expect_false(is.na(out$pi_power))
})

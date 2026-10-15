# The SIM_* functions forward 'nonborrow_prior' to their SCI_* fit. The default
# must reproduce the previous behaviour, and a custom prior must change the
# operating characteristics. EB path + parallel = FALSE keep these fast (no
# gMAP, no cluster). The normal endpoint uses the same one-line pass-through and
# is covered by the SCI_norm helper tests in test-nonborrow-prior.R.

test_that("SIM_pois forwards nonborrow_prior (default = old behaviour, custom bites)", {
        a <- list(n_sim = 6, H = 6, offset = 3, lambda = 8, shape = 12,
                  shift = c(0, 0), prior_est = "EB", seed = 1, parallel = FALSE)
        base0 <- suppressWarnings(do.call(SIM_pois, a))
        base1 <- suppressWarnings(do.call(SIM_pois, c(a, list(nonborrow_prior = c(0.5, 0.001)))))
        cust  <- suppressWarnings(do.call(SIM_pois, c(a, list(nonborrow_prior = c(5, 1)))))
        expect_equal(base0, base1)
        expect_false(isTRUE(all.equal(base0, cust)))
})

test_that("SIM_binom forwards nonborrow_prior (default = Beta(1,1), custom bites)", {
        a <- list(n_sim = 6, H = 5, n = 40, pi = 0.3, rho = 0.05,
                  shift = c(0, 0), prior_est = "EB", seed = 1, parallel = FALSE)
        base0 <- suppressWarnings(do.call(SIM_binom, a))
        base1 <- suppressWarnings(do.call(SIM_binom, c(a, list(nonborrow_prior = c(1, 1)))))
        cust  <- suppressWarnings(do.call(SIM_binom, c(a, list(nonborrow_prior = c(0.5, 0.5)))))
        expect_equal(base0, base1)
        expect_false(isTRUE(all.equal(base0, cust)))
})

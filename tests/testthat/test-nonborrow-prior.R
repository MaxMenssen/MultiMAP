# The 'nonborrow_prior' argument exposes the previously hardcoded prior for the
# non-borrowing arms. The default must reproduce the old behaviour exactly, a
# custom prior must change the non-borrowing posteriors, and bad input must be
# rejected. These use the EB path / the pure-arithmetic normal helper, so they
# run fast (no gMAP).

test_that("SCI_pois: default nonborrow_prior reproduces Gamma(0.5, 0.001)", {
        hp <- data.frame(events = c(2, 1, 3, 2, 1, 2), offset = 5)
        np <- data.frame(group = c("0", "10", "50"), events = c(2, 1, 4), offset = 5)
        a <- suppressWarnings(SCI_pois(hp, np, prior_est = "EB", seed = 1, n_mcmc = 3000))
        b <- suppressWarnings(SCI_pois(hp, np, prior_est = "EB", seed = 1, n_mcmc = 3000,
                                       nonborrow_prior = c(0.5, 0.001)))
        expect_equal(a$comp, b$comp)

        d <- suppressWarnings(SCI_pois(hp, np, prior_est = "EB", seed = 1, n_mcmc = 3000,
                                       nonborrow_prior = c(3, 1)))
        expect_false(isTRUE(all.equal(a$comp$estimate, d$comp$estimate)))

        expect_error(SCI_pois(hp, np, nonborrow_prior = c(-1, 1)))
        expect_error(SCI_pois(hp, np, nonborrow_prior = 0.5))
})

test_that("SCI_binom: default nonborrow_prior reproduces Beta(1, 1)", {
        hb <- data.frame(ev = c(5, 3, 4, 6, 2), nonev = c(20, 22, 21, 19, 23))
        nb <- data.frame(group = c("0", "10", "50"), ev = c(4, 6, 9), nonev = c(21, 19, 16))
        a <- suppressWarnings(SCI_binom(hb, nb, prior_est = "EB", seed = 1, n_mcmc = 3000))
        b <- suppressWarnings(SCI_binom(hb, nb, prior_est = "EB", seed = 1, n_mcmc = 3000,
                                        nonborrow_prior = c(1, 1)))
        expect_equal(a$comp, b$comp)

        d <- suppressWarnings(SCI_binom(hb, nb, prior_est = "EB", seed = 1, n_mcmc = 3000,
                                        nonborrow_prior = c(0.5, 0.5)))
        expect_false(isTRUE(all.equal(a$comp$estimate, d$comp$estimate)))

        expect_error(SCI_binom(hb, nb, nonborrow_prior = c(0, 1)))
        expect_error(SCI_binom(hb, nb, nonborrow_prior = c(1, 1, 1)))
})

test_that("SCI_norm: NULL nonborrow_prior is the flat-prior posterior, c(m, s) is conjugate", {
        # NULL -> improper flat prior -> posterior N(ybar, se)
        flat <- .norm_nonborrow_post(ybar = c(10, 12), se = c(0.5, 0.4), nonborrow_prior = NULL)
        expect_equal(flat$m, c(10, 12))
        expect_equal(flat$s, c(0.5, 0.4))

        # proper N(m0, s0): posterior mean between m0 and ybar, posterior sd < se
        prp <- .norm_nonborrow_post(ybar = 12, se = 0.5, nonborrow_prior = c(0, 1))
        expect_true(prp$m > 0 && prp$m < 12)
        expect_lt(prp$s, 0.5)
        # closed-form check
        prec <- 1 / 1^2 + 1 / 0.5^2
        expect_equal(prp$m, (0 / 1 + 12 / 0.25) / prec)
        expect_equal(prp$s, sqrt(1 / prec))

        hn <- data.frame(mean = c(10, 11, 9), sd = c(2, 2.2, 1.9), n = c(30, 28, 32))
        nn <- data.frame(group = c("0", "10"), mean = c(10, 12), sd = c(2, 2), n = c(30, 30))
        expect_error(SCI_norm(hn, nn, nonborrow_prior = c(0, -1)))
        expect_error(SCI_norm(hn, nn, nonborrow_prior = 5))
})

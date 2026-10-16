# The summary() ESS block previously always said "ESS MAP", even when borrowing
# used the empirical-Bayes prior (prior_est = "EB"). The object now stores
# prior_est and summary() labels the effective sample sizes with the method
# actually used. EB path -> fast, no gMAP; the MAP check needs gMAP.

test_that("summary() labels EB borrowing as EB, not MAP", {
        hb <- data.frame(events = c(2, 2, 1, 5, 3, 5), non_events = c(48, 47, 49, 44, 51, 44))
        nb <- data.frame(group = c("0", "10", "50"), events = c(3, 4, 8), non_events = c(44, 36, 39))

        reb <- suppressWarnings(SCI_binom(hb, nb, prior_est = "EB", n_mcmc = 1500, seed = 1))
        expect_identical(reb$prior_est, "EB")
        out <- capture.output(summary(reb))
        expect_true(any(grepl("ESS  EB:", out, fixed = TRUE)))
        expect_false(any(grepl("ESS  MAP:", out, fixed = TRUE)))

        # PI_binom and the multi-arm path too
        peb <- suppressWarnings(PI_binom(hb, n0 = 50, prior_est = "EB", n_mcmc = 1500, seed = 1))
        expect_true(any(grepl("ESS  EB:", capture.output(summary(peb)), fixed = TRUE)))

        meb <- suppressWarnings(SCI_binom(list("0" = hb), nb, prior_est = "EB", n_mcmc = 1500, seed = 1))
        expect_true(any(grepl("ESS  EB:", capture.output(summary(meb)), fixed = TRUE)))

        # Poisson EB as well
        hp <- data.frame(events = c(10, 12, 8, 15, 11), offset = 5)
        np <- data.frame(group = c("0", "100"), events = c(10, 15), offset = 5)
        rp <- suppressWarnings(SCI_pois(hp, np, prior_est = "EB", n_mcmc = 1500, seed = 1))
        expect_true(any(grepl("ESS  EB:", capture.output(summary(rp)), fixed = TRUE)))
})

test_that("summary() labels MAP borrowing as MAP", {
        skip_on_cran()
        hb <- data.frame(events = c(2, 2, 1, 5, 3, 5), non_events = c(48, 47, 49, 44, 51, 44))
        nb <- data.frame(group = c("0", "10", "50"), events = c(3, 4, 8), non_events = c(44, 36, 39))
        rmap <- suppressWarnings(SCI_binom(hb, nb, prior_est = "MAP", n_mcmc = 1500, seed = 1))
        expect_identical(rmap$prior_est, "MAP")
        out <- capture.output(summary(rmap))
        expect_true(any(grepl("ESS  MAP:", out, fixed = TRUE)))
        expect_false(any(grepl("ESS  EB:", out, fixed = TRUE)))
})

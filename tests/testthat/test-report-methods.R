# The REPORT_*() reports now include a Methods section that depends on the input:
# the prior method (MAP vs EB), its hyperpriors (MAP) or estimated prior and
# overdispersion (EB), and the robustification. This needs the EB overdispersion
# rho to survive in the object (eb_estimates). EB path -> fast; MAP needs gMAP.

test_that("eb_estimates preserves the EB overdispersion rho (binomial)", {
        hb <- data.frame(events = c(2, 2, 1, 5, 3, 5), non_events = c(48, 47, 49, 44, 51, 44))
        nb <- data.frame(group = c("0", "10"), events = c(3, 4), non_events = c(44, 36))

        eb <- suppressWarnings(SCI_binom(hb, nb, prior_est = "EB", n_mcmc = 1500, seed = 1))
        expect_true(all(c("pi", "rho") %in% names(eb$eb_estimates)))
        expect_true(is.finite(eb$eb_estimates[["rho"]]))

        pe <- suppressWarnings(PI_binom(hb, n0 = 50, prior_est = "EB", n_mcmc = 1500, seed = 1))
        expect_true("rho" %in% names(pe$eb_estimates))

        me <- suppressWarnings(SCI_binom(list("0" = hb), nb, prior_est = "EB", n_mcmc = 1500, seed = 1))
        expect_true("rho" %in% names(me$eb_estimates[["0"]]))
})

test_that(".report_methods() reflects the prior method and endpoint", {
        hb <- data.frame(events = c(2, 2, 1, 5, 3, 5), non_events = c(48, 47, 49, 44, 51, 44))
        nb <- data.frame(group = c("0", "10"), events = c(3, 4), non_events = c(44, 36))
        eb <- suppressWarnings(SCI_binom(hb, nb, prior_est = "EB", n_mcmc = 1500, seed = 1))
        m  <- .report_methods(eb)
        expect_true(grepl("empirical-Bayes", m))
        expect_true(grepl("overdispersion rho", m))      # rho reported
        expect_true(grepl("Lui, 2000", m))
        expect_false(grepl("meta-analytic", m))

        hp <- data.frame(events = c(10, 12, 8, 15, 9, 11), offset = 4)
        np <- data.frame(group = c("0", "10"), events = c(11, 14), offset = 4)
        mp <- .report_methods(suppressWarnings(SCI_pois(hp, np, prior_est = "EB", n_mcmc = 1500, seed = 1)))
        expect_true(grepl("Gamma\\(", mp))
        expect_true(grepl("Tarone, 1982", mp))
})

test_that("MAP methods paragraph names the hyperpriors", {
        skip_on_cran()
        hb <- data.frame(events = c(2, 2, 1, 5, 3, 5), non_events = c(48, 47, 49, 44, 51, 44))
        nb <- data.frame(group = c("0", "10"), events = c(3, 4), non_events = c(44, 36))
        mp <- suppressWarnings(SCI_binom(hb, nb, prior_est = "MAP", n_mcmc = 1500, seed = 1))
        expect_true(all(is.na(mp$eb_estimates)))
        m <- .report_methods(mp)
        expect_true(grepl("meta-analytic-predictive", m))
        expect_true(grepl("HalfNormal\\(0,", m))         # between-study heterogeneity prior
        expect_true(grepl("prior on the intercept", m))
        expect_false(grepl("empirical-Bayes", m))
})

test_that("REPORT_binom embeds a Methods section reflecting prior_est", {
        hb <- data.frame(events = c(2, 2, 1, 5, 3, 5), non_events = c(48, 47, 49, 44, 51, 44))
        nb <- data.frame(group = c("0", "10"), events = c(3, 4), non_events = c(44, 36))
        eb <- suppressWarnings(SCI_binom(hb, nb, prior_est = "EB", n_mcmc = 1500, seed = 1))

        qmd <- tempfile(fileext = ".qmd")
        suppressMessages(REPORT_binom(eb, file = qmd))
        on.exit(unlink(c(qmd, sub("\\.qmd$", "_object.rds", qmd))), add = TRUE)

        L <- readLines(qmd)
        expect_true(any(grepl("^# Methods", L)))
        expect_true(any(grepl("empirical-Bayes", L)))
        expect_true(any(grepl("overdispersion rho", L)))
})

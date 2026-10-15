# Regression (issue H4): data-frame inputs are documented by column name but
# were read by position (and the names overwritten by position), so correctly
# named columns supplied in a different order gave silently wrong results.
# Columns are now selected by name when all expected names are present, with a
# positional fallback (and a message) otherwise. EB paths need no gMAP; the
# normal endpoint is MAP-only, so that test skips on CRAN.

test_that("SCI_binom / PI_binom read histdat columns by name, not position", {
        hb <- data.frame(events = c(2, 3, 1, 4, 2), non_events = c(48, 47, 49, 46, 48))
        nb <- data.frame(group = c("control", "T1", "T2"),
                         events = c(2, 5, 9), non_events = c(48, 45, 41))

        a <- suppressWarnings(SCI_binom(hb, nb, "EB", n_mcmc = 2000, seed = 1))
        b <- suppressWarnings(SCI_binom(hb[, c("non_events", "events")], nb,
                                        "EB", n_mcmc = 2000, seed = 1))
        expect_equal(a$pred_int, b$pred_int, tolerance = 1e-6)
        expect_equal(a$comp,     b$comp,     tolerance = 1e-6)

        p  <- suppressWarnings(PI_binom(hb, n0 = 50, prior_est = "EB", n_mcmc = 2000, seed = 1))
        p2 <- suppressWarnings(PI_binom(hb[, c("non_events", "events")], n0 = 50,
                                        prior_est = "EB", n_mcmc = 2000, seed = 1))
        expect_equal(p$pred_int, p2$pred_int, tolerance = 1e-6)
})

test_that("SCI_pois reads newdat columns by name, not position", {
        hp <- data.frame(events = c(10, 12, 8, 15, 11), offset = c(50, 55, 45, 60, 50))
        np <- data.frame(group = c("control", "T1", "T2"),
                         events = c(10, 15, 20), offset = c(50, 52, 48))

        a <- suppressWarnings(SCI_pois(hp, np, "EB", n_mcmc = 2000, seed = 1))
        b <- suppressWarnings(SCI_pois(hp, np[, c("group", "offset", "events")],
                                       "EB", n_mcmc = 2000, seed = 1))
        expect_equal(a$comp, b$comp, tolerance = 1e-6)
})

test_that("SCI_norm reads columns by name and labels $newdat correctly", {
        skip_on_cran()
        hn <- data.frame(mean = c(9.8, 10.2, 9.5, 10, 10.4, 9.9),
                         sd = c(2.1, 1.9, 2.3, 2, 2.2, 2), n = c(40, 55, 38, 47, 20, 30))
        nn <- data.frame(group = c("control", "d1", "d2"),
                         mean = c(10.1, 11.4, 12.7), sd = c(2, 2.2, 2.1), n = c(30, 30, 30))

        a <- suppressWarnings(SCI_norm(hn, nn, n_mcmc = 2000, seed = 1))
        b <- suppressWarnings(SCI_norm(hn, nn[, c("group", "n", "mean", "sd")],
                                       n_mcmc = 2000, seed = 1))
        expect_equal(a$comp, b$comp, tolerance = 1e-6)
        # $newdat must carry the true means, not the reordered 'n' column
        expect_equal(b$newdat$mean, nn$mean)
})

test_that("positional fallback is used (with a message) only when names are absent", {
        # non-canonical newdat names -> positional fallback + message (histdat is
        # canonical here, so exactly one message fires)
        hb0 <- data.frame(events = c(2, 3, 1, 4, 2), non_events = c(48, 47, 49, 46, 48))
        nx  <- data.frame(g = c("control", "T1"), e = c(2, 5), ne = c(48, 45))
        expect_message(suppressWarnings(SCI_binom(hb0, nx, "EB", n_mcmc = 2000, seed = 1)),
                       "by position")

        # canonical names -> no positional-fallback message
        hb <- data.frame(events = c(2, 3, 1, 4, 2), non_events = c(48, 47, 49, 46, 48))
        nb <- data.frame(group = c("control", "T1"), events = c(2, 5), non_events = c(48, 45))
        expect_no_message(suppressWarnings(SCI_binom(hb, nb, "EB", n_mcmc = 2000, seed = 1)))
})

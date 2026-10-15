# Regression (issue M6): SCI_* / PI_* validated only the shape of histdat /
# newdat, never the values. NA / negative / non-integer counts, zero or negative
# offsets, and non-positive sd / n were silently used (producing all-NA or
# degenerate SCIs) or failed with cryptic internal errors. Values are now
# validated at entry with explicit messages. EB path -> fast, no gMAP.

test_that("SCI_binom rejects invalid count values (histdat and newdat)", {
        hb <- data.frame(events = c(2, 2, 1, 5), non_events = c(48, 47, 49, 44))
        nb <- data.frame(group = c("0", "10"), events = c(3, 4), non_events = c(44, 36))

        bad <- nb; bad$events[2] <- NA
        expect_error(SCI_binom(hb, bad, prior_est = "EB", n_mcmc = 1000, seed = 1),
                     "missing or non-finite")
        bad <- nb; bad$events[2] <- -1
        expect_error(SCI_binom(hb, bad, prior_est = "EB", n_mcmc = 1000, seed = 1),
                     "non-negative")
        bad <- nb; bad$events[2] <- 4.5
        expect_error(SCI_binom(hb, bad, prior_est = "EB", n_mcmc = 1000, seed = 1),
                     "whole numbers")
        bad <- hb; bad$non_events[1] <- NA
        expect_error(SCI_binom(bad, nb, prior_est = "EB", n_mcmc = 1000, seed = 1),
                     "missing or non-finite")
        # valid data is unaffected
        expect_error(suppressWarnings(SCI_binom(hb, nb, prior_est = "EB", n_mcmc = 1000, seed = 1)), NA)
})

test_that("SCI_pois rejects non-positive offsets, PI_pois too", {
        hp <- data.frame(events = c(10, 12, 8, 15), offset = c(50, 55, 45, 60))
        np <- data.frame(group = c("0", "a"), events = c(11, 14), offset = 50)

        bad <- hp; bad$offset[1] <- 0
        expect_error(SCI_pois(bad, np, prior_est = "EB", n_mcmc = 1000, seed = 1), "positive")
        bad <- np; bad$events[1] <- -1
        expect_error(SCI_pois(hp, bad, prior_est = "EB", n_mcmc = 1000, seed = 1), "non-negative")
        bad <- hp; bad$offset[2] <- 0
        expect_error(PI_pois(bad, offset0 = 50, prior_est = "EB", n_mcmc = 1000, seed = 1), "positive")
})

test_that("SCI_norm rejects non-positive sd / n, PI_norm too", {
        hn <- data.frame(mean = c(1.1, 1, 1.37, 1.19), sd = c(.32, .3, .37, .33), n = c(40, 5, 5, 5))
        nn <- data.frame(group = c("0", "a"), mean = c(1.05, 1.84), sd = c(.31, .46), n = 5)

        bad <- hn; bad$n[1] <- -40
        expect_error(SCI_norm(bad, nn, n_mcmc = 1000, seed = 1), "positive")
        bad <- nn; bad$sd[1] <- 0
        expect_error(SCI_norm(hn, bad, n_mcmc = 1000, seed = 1), "positive")
        bad <- hn; bad$sd[2] <- 0
        expect_error(PI_norm(bad, n0 = 30, seed = 1), "positive")
})

test_that("multi-arm paths validate each arm's historical data", {
        hb <- data.frame(events = c(2, 2, 1, 5), non_events = c(48, 47, 49, 44))
        bad <- data.frame(events = c(2, -1, 1, 5), non_events = c(48, 47, 49, 44))
        nb <- data.frame(group = c("0", "10"), events = c(3, 4), non_events = c(44, 36))
        expect_error(SCI_binom(list("0" = bad), nb, prior_est = "EB", n_mcmc = 1000, seed = 1),
                     "non-negative")
})

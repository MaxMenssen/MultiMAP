# Regression (issue H3): the binomial empirical-Bayes "not overdispersed" cap
# compared b_est with sum(events) instead of sum(non-events). For a pooled
# proportion above 0.5 the fitted Beta prior could then carry several times more
# information (a + b) than all the historical data (here a + b = 918.7 from
# N = 500), over-shrinking the concurrent control, and the result depended on
# which outcome was coded as the event. The cap now triggers on total
# information (a + b) > N and falls back to the pooled Beta(sum ev, sum non-ev).
# EB path only (no gMAP), so these run fast.

test_that("binomial EB prior never carries more information than the historical data", {
        hd <- data.frame(events = c(44, 44, 46, 47, 39, 46, 46, 45, 45, 43))
        hd$non_events <- 50 - hd$events
        N  <- sum(hd$events) + sum(hd$non_events)          # 500
        ab <- function(m) sum(m[c("a", "b"), 1])           # prior pseudo-count a + b

        # multi-arm helper .map_prior_binom(), was Beta(817.7, 101.1) -> a+b 918.7
        expect_lte(ab(suppressWarnings(
                .map_prior_binom(hd, "EB", NULL, "HalfNormal", NULL)$map)), N)

        # single-arm PI_binom() inline EB path
        expect_lte(ab(suppressWarnings(
                PI_binom(hd, n0 = 50, prior_est = "EB", n_mcmc = 2000, seed = 1)$prior)), N)

        # single-arm SCI_binom() inline EB path
        nd <- data.frame(group = c("C", "T"), events = c(45, 40), non_events = c(5, 10))
        expect_lte(ab(suppressWarnings(
                SCI_binom(hd, nd, prior_est = "EB", contr = "pi_diff",
                          n_mcmc = 2000, seed = 1)$prior)), N)
})

test_that("binomial EB prior does not depend on event/non-event coding", {
        hd <- data.frame(events = c(44, 44, 46, 47, 39, 46, 46, 45, 45, 43))
        hd$non_events <- 50 - hd$events
        hd2 <- data.frame(events = hd$non_events, non_events = hd$events)

        m  <- suppressWarnings(.map_prior_binom(hd,  "EB", NULL, "HalfNormal", NULL)$map)
        m2 <- suppressWarnings(.map_prior_binom(hd2, "EB", NULL, "HalfNormal", NULL)$map)
        # recoding events <-> non-events must simply swap (a, b)
        expect_equal(unname(m[c("a", "b"), 1]), unname(m2[c("b", "a"), 1]), tolerance = 1e-6)
})

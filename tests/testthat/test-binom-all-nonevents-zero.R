# Regression (issue M1): historical binomial data with no non-events (a 100%
# background rate) only had a special case for the mirror situation (all events
# zero). In the EB path the moment estimator then produced the degenerate
# Beta(sum events, 0), which (after the H2 ESS guard) ran but gave a point mass
# at rate 1 and an NA effective sample size. The EB path now mirrors the
# all-events-zero branch: Beta(N - 1, 1) with N = total observations. The MAP
# path is already handled (fixed tau = 1, issue H1).

test_that("binomial EB: all non-events zero gives Beta(N-1, 1), mirroring all-events-zero", {
        hd <- data.frame(events = c(50, 40, 45), non_events = c(0, 0, 0))   # 100% rate
        N  <- sum(hd$events) + sum(hd$non_events)                           # 135
        ab <- function(m) unname(m[c("a", "b"), 1])

        # multi-arm helper .map_prior_binom(): was the degenerate Beta(135, 0)
        expect_equal(ab(suppressWarnings(
                .map_prior_binom(hd, "EB", NULL, "HalfNormal", NULL)$map)), c(N - 1, 1))

        # single-arm SCI_binom(): valid prior (both params > 0) and finite pred_int
        nd <- data.frame(group = c("C", "T"), events = c(45, 40), non_events = c(5, 10))
        r  <- suppressWarnings(SCI_binom(hd, nd, prior_est = "EB", n_mcmc = 2000, seed = 1))
        expect_equal(ab(r$prior), c(N - 1, 1))
        expect_true(all(is.finite(r$pred_int)))

        # single-arm PI_binom()
        p  <- suppressWarnings(PI_binom(hd, n0 = 50, prior_est = "EB", n_mcmc = 2000, seed = 1))
        expect_equal(ab(p$prior), c(N - 1, 1))
        expect_true(all(is.finite(p$pred_int)))

        # the mirror case is unchanged: all events zero still gives Beta(1, N-1)
        hz <- data.frame(events = c(0, 0, 0), non_events = c(50, 40, 45))
        expect_equal(ab(suppressWarnings(
                .map_prior_binom(hz, "EB", NULL, "HalfNormal", NULL)$map)), c(1, N - 1))
})

test_that("binomial MAP: all non-events zero runs and gives a prior near rate 1", {
        skip_on_cran()
        hd <- data.frame(events = c(50, 40, 45), non_events = c(0, 0, 0))
        nd <- data.frame(group = c("C", "T"), events = c(45, 40), non_events = c(5, 10))
        r  <- suppressWarnings(SCI_binom(hd, nd, prior_est = "MAP", n_mcmc = 2000, seed = 1))
        expect_true(all(is.finite(r$pred_int)))
        expect_gt(as.numeric(summary(r$prior)["mean"]), 0.9)
})

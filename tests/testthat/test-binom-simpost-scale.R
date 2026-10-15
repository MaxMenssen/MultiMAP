# Regression (issue M2): SCI_binom single-arm with contr = "OR" overwrote the
# posterior-sample matrix with its logit and returned that as $simpost, while
# contr = "RR"/"pi_diff" and the whole multi-arm path keep $simpost on the
# proportion scale. $comp was unaffected, but user post-processing of $simpost
# silently got logit-scale values in this one configuration. The column names
# also differed. The logit transform is now applied locally, and the single-arm
# columns are named by group like the multi-arm path. EB path -> fast, no gMAP.

test_that("SCI_binom OR keeps $simpost on the proportion scale with group-named columns", {
        hb <- data.frame(events = c(2, 2, 1, 5, 3, 5, 2, 2),
                         non_events = c(48, 47, 49, 44, 51, 44, 54, 53))
        nb <- data.frame(group = c("0", "10", "50"),
                         events = c(3, 4, 8), non_events = c(44, 36, 39))

        s <- suppressWarnings(SCI_binom(hb, nb, prior_est = "EB", contr = "OR",
                                        seed = 1, n_mcmc = 2000))
        # proportions in [0, 1], not logits (which are negative for rare events)
        expect_gte(min(s$simpost), 0)
        expect_lte(max(s$simpost), 1)
        # columns named by group, consistent with the multi-arm path
        expect_equal(colnames(s$simpost), as.character(nb$group))

        # the contrast scale is unchanged: $simpost no longer depends on 'contr'
        r <- suppressWarnings(SCI_binom(hb, nb, prior_est = "EB", contr = "RR",
                                        seed = 1, n_mcmc = 2000))
        expect_equal(s$simpost, r$simpost)
})

test_that("SCI_binom multi-arm OR $simpost is proportion-scale with group-named columns", {
        hb <- data.frame(events = c(2, 2, 1, 5, 3, 5, 2, 2),
                         non_events = c(48, 47, 49, 44, 51, 44, 54, 53))
        nb <- data.frame(group = c("0", "10", "50"),
                         events = c(3, 4, 8), non_events = c(44, 36, 39))

        m <- suppressWarnings(SCI_binom(list("0" = hb), nb, prior_est = "EB",
                                        contr = "OR", seed = 1, n_mcmc = 2000))
        expect_gte(min(m$simpost), 0)
        expect_lte(max(m$simpost), 1)
        expect_equal(colnames(m$simpost), as.character(nb$group))
})

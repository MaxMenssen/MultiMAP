# Regression (issue M7): the binomial MAP default intercept prior was N(0, 2) on
# the logit scale, centred at p = 0.5, unlike the data-centred normal (historical
# mean) and Poisson (log lambda_bar) families. For rare / frequent events it
# pulled the MAP prior towards 0.5. The default is now Normal(logit(pooled rate),
# 2) with a continuity correction. Needs gMAP, so skip on CRAN.

test_that("binomial MAP default intercept is data-centred, not logit-0", {
        skip_on_cran()
        hd <- data.frame(events = c(0, 1, 0), non_events = 50 - c(0, 1, 0))   # pooled ~0.0067
        pm <- function(m) sum(m["w", ] * m["a", ] / (m["a", ] + m["b", ]))
        p_pool_cc <- (1 + 0.5) / (150 + 1)

        # multi-arm helper .map_prior_binom(): the NULL default now matches an
        # explicit data-centred intercept (before the fix it matched logit 0)
        set.seed(11); d_def <- pm(suppressWarnings(
                .map_prior_binom(hd, "MAP", NULL, "HalfNormal", NULL))$map)
        set.seed(11); d_ctr <- pm(suppressWarnings(
                .map_prior_binom(hd, "MAP", NULL, "HalfNormal", cbind(qlogis(p_pool_cc), 2)))$map)
        expect_equal(d_def, d_ctr, tolerance = 1e-6)
        expect_lt(d_def, 0.02)                 # close to the pooled rate, not pulled to ~0.027

        # single-arm SCI_binom() and PI_binom() use the same data-centred default
        s <- suppressWarnings(SCI_binom(hd, data.frame(group = c("0", "1"),
                                                       events = c(0, 2), non_events = c(50, 48)),
                                        prior_est = "MAP", seed = 11))
        expect_lt(pm(s$prior), 0.02)
        p <- suppressWarnings(PI_binom(hd, n0 = 50, prior_est = "MAP", seed = 11))
        expect_lt(pm(p$prior), 0.02)
})

test_that("all-success historical data does not crash the MAP intercept (continuity correction)", {
        skip_on_cran()
        hs <- data.frame(events = c(50, 49, 50), non_events = c(0, 1, 0))   # ~100% rate
        pm <- function(m) sum(m["w", ] * m["a", ] / (m["a", ] + m["b", ]))
        set.seed(11)
        m <- suppressWarnings(.map_prior_binom(hs, "MAP", NULL, "HalfNormal", NULL)$map)
        expect_true(is.finite(pm(m)))
        expect_gt(pm(m), 0.9)                  # centred near 1, not pulled down to 0.5
})

# Regression (issue H2): a numerically unstable MAP fit (a degenerate automixfit
# mixture, e.g. a component with an overflowed parameter from a very diffuse
# tau.prior) makes BOTH ess() estimators fail. The elir estimator was guarded by
# try(), but the method-of-moments fallback was not, so PI_binom / SCI_binom (and
# the pois / norm analogues) aborted with the cryptic "missing value where
# TRUE/FALSE needed" / "assert_that: missing values present". .ess_safe() now
# guards both and returns NA with an informative warning. Fast: no gMAP.

test_that(".ess_safe returns NA with an informative warning on a degenerate mixture", {
        # Beta component with an overflowed shape -> both elir and moment fail.
        degen <- RBesT::mixbeta(c(1, 1e308, 1))

        expect_warning(res <- .ess_safe(degen, cap = 100), "numerically unstable")
        expect_true(is.na(res))
})

test_that(".ess_safe still returns a finite, capped ESS for a healthy mixture", {
        m <- RBesT::mixbeta(c(1, 20, 80))   # Beta(20, 80): ESS ~ a + b = 100

        res <- expect_no_warning(.ess_safe(m, cap = 100))
        expect_true(is.finite(res))
        expect_gte(res, 0)
        expect_lte(res, 100)
})

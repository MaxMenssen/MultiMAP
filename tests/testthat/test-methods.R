# Shared S3 methods, exercised on a fast EB binomial fit.

fit_eb <- function(contr = "pi_diff", limits = "two.sided"){
        histdat <- data.frame(events     = c(11, 8, 14, 9),
                              non_events = c(29, 32, 26, 31))
        newdat  <- data.frame(group      = c("control", "d1", "d454"),
                              events     = c(8, 14, 17),
                              non_events = c(22, 16, 13))
        suppressWarnings(SCI_binom(histdat, newdat, prior_est = "EB", contr = contr,
                                   limits = limits, seed = 1))
}

test_that("plot.MultiMAP returns ggplot objects for all three views (binomial)", {
        res <- fit_eb()
        expect_s3_class(plot(res, which = 1), "ggplot")
        expect_s3_class(plot(res, which = 2), "ggplot")
        expect_s3_class(plot(res, which = 3), "ggplot")
})

test_that("plot which = 2 title and null value follow the contrast type", {
        expect_s3_class(plot(fit_eb(contr = "RR"), which = 2), "ggplot")
        expect_s3_class(plot(fit_eb(contr = "OR"),      which = 2), "ggplot")
})

test_that("summary and print run and return invisibly (binomial)", {
        res <- fit_eb()
        expect_output(summary(res), "proportion")
        expect_output(summary(res), "ratio of proportions|difference of proportions")
        expect_output(print(res), "control")
        expect_invisible(summary(res))
        expect_invisible(print(res))
})

test_that(".null_value maps ratio-type contrasts to 1 and differences to 0", {
        nv <- getFromNamespace(".null_value", "MultiMAP")
        expect_equal(nv("pi_diff"),    0)
        expect_equal(nv("RR"),    1)
        expect_equal(nv("OR"),         1)
        expect_equal(nv("mean_diff"),  0)
        expect_equal(nv("mean_ratio"), 1)
})

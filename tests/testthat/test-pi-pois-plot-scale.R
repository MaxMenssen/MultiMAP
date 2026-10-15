# Regression (issue H8): plot() of a PI_pois object drew the gamma-Poisson
# predictive on the raw COUNT scale (because offset0 was not stored, plot used
# exposure = 1) while labelling the x-axis "rate" and comparing the counts with
# the rate-scale pred_int, so the highlighted "inside PI" bars sat in the far
# left tail with almost no mass. PI_pois now stores offset0 and the plot uses it
# as the exposure. EB path -> no gMAP, fast.

test_that("plot(PI_pois) is drawn on the rate scale and highlights the PI", {
        hp <- data.frame(events = c(18, 22, 15, 25, 20), offset = 5)
        pp <- suppressWarnings(PI_pois(hp, offset0 = 5, prior_est = "EB", seed = 1))

        expect_equal(unname(pp$offset0), 5)                 # offset0 is stored

        b  <- ggplot2::ggplot_build(plot(pp))$data[[1]]
        hl <- b$fill == "lightblue"

        # the mode is at the rate (pred_int median ~ 4), not the raw count (~ 20)
        expect_lt(b$x[which.max(b$y)], 10)
        # highlighted bars coincide with pred_int [lower, upper] on the rate scale
        expect_equal(min(b$x[hl]), unname(pp$pred_int["lower"]), tolerance = 0.3)
        expect_equal(max(b$x[hl]), unname(pp$pred_int["upper"]), tolerance = 0.3)
        # and they carry most of the probability mass, not a far tail
        expect_gt(sum(b$y[hl]), 0.5)
})

test_that("plot(PI_pois) multi-arm uses the per-arm offset0 in each facet", {
        hp <- data.frame(events = c(18, 22, 15, 25, 20), offset = 5)
        pm <- suppressWarnings(PI_pois(
                list(A = hp, B = data.frame(events = c(40, 44, 38, 50), offset = 10)),
                offset0 = c(A = 5, B = 10), prior_est = "EB", seed = 1))

        expect_equal(pm$offset0[["A"]], 5)
        expect_equal(pm$offset0[["B"]], 10)

        bd <- ggplot2::ggplot_build(plot(pm))$data[[1]]
        # every panel's mode is on the rate scale (< 10), not the count scale (~20/42)
        modes <- tapply(seq_len(nrow(bd)), bd$PANEL,
                        function(i) bd$x[i][which.max(bd$y[i])])
        expect_true(all(modes < 10))
})

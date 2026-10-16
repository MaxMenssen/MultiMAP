# Single-arm borrowing: plot(res, which = 3) previously showed the prior /
# likelihood / posterior only for the borrowing (control) arm; the non-borrowing
# treatment arms were missing. It now facets over all arms, like the multi-arm
# which = 3 plot, with the treatment arms marked "(no borrowing)". EB path ->
# fast, no gMAP.

n_panels <- function(p) length(unique(ggplot2::ggplot_build(p)$data[[1]]$PANEL))

test_that("single-arm SCI_binom which = 3 facets over all arms (incl. non-borrowing)", {
        hb <- data.frame(events = c(2, 2, 1, 5, 3, 5), non_events = c(48, 47, 49, 44, 51, 44))
        nb <- data.frame(group = c("0", "10", "50"), events = c(3, 4, 8), non_events = c(44, 36, 39))
        r  <- suppressWarnings(SCI_binom(hb, nb, prior_est = "EB", n_mcmc = 1500, seed = 1))

        p  <- plot(r, which = 3)
        expect_s3_class(p, "ggplot")
        expect_equal(n_panels(p), nrow(nb))                 # one panel per arm, not just the control

        lv <- levels(ggplot2::ggplot_build(p)$plot$data$arm)
        # the two treatment arms are shown and flagged as non-borrowing
        expect_equal(sum(grepl("no borrowing", lv)), nrow(nb) - 1L)
        expect_true("0" %in% lv)                            # the borrowing control keeps its plain label

        # which = 1 / which = 2 are unchanged (single panel)
        expect_equal(n_panels(plot(r, which = 1)), 1L)
        expect_equal(n_panels(plot(r, which = 2)), 1L)
})

test_that("single-arm SCI_pois which = 3 also facets over all arms", {
        hp <- data.frame(events = c(10, 12, 8, 15, 11), offset = 5)
        np <- data.frame(group = c("0", "100", "200"), events = c(10, 15, 20), offset = 5)
        r  <- suppressWarnings(SCI_pois(hp, np, prior_est = "EB", n_mcmc = 1500, seed = 1))
        expect_equal(n_panels(plot(r, which = 3)), nrow(np))
})

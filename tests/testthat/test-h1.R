test_that("h1_norm returns the expected structure and arms", {
        set.seed(1)
        d <- h1_norm(mu = 10, H = 5, n = 30, sd_h = 2, sd_hi = 3,
                     shift = c(0, 1.5))
        expect_s3_class(d, "data.frame")
        expect_named(d, c("y_hi", "a", "arm", "e"))
        # 5 control blocks x 3 arms x 30 obs
        expect_equal(nrow(d), 5 * 3 * 30)
        expect_setequal(levels(d$arm), c("control", "trt1", "trt2"))
        expect_equal(nlevels(d$a), 5)
})

test_that("h1_norm shift = NULL gives controls only", {
        set.seed(1)
        d <- h1_norm(mu = 10, H = 3, n = 20, sd_h = 2, sd_hi = 3, shift = NULL)
        expect_setequal(levels(d$arm), "control")
        expect_equal(nrow(d), 3 * 20)
})

test_that("h1_binom returns per-study/arm counts with shared cluster proportion", {
        set.seed(1)
        d <- h1_binom(pi = 0.3, H = 4, n = 40, rho = 0.05)
        expect_s3_class(d, "data.frame")
        expect_named(d, c("study", "arm", "events", "n", "pi_cluster", "prob"))
        expect_equal(nrow(d), 4)                    # control only, 4 studies
        expect_true(all(d$arm == "control"))
        expect_true(all(d$n == 40))
        expect_true(all(d$events >= 0 & d$events <= 40))
})

test_that("h1_binom treatment arms share the cluster proportion", {
        set.seed(2)
        d <- h1_binom(pi = 0.3, H = 1, n = 40, rho = 0.05, shift = c(0, 0.15))
        expect_setequal(as.character(d$arm), c("control", "trt1", "trt2"))
        # control and trt1 (shift 0) share pi_cluster; on prob scale prob == pi_cluster
        expect_equal(d$prob[d$arm == "control"], d$pi_cluster[d$arm == "control"])
        expect_equal(d$prob[d$arm == "trt1"],    d$pi_cluster[d$arm == "trt1"])
        # trt2 offset by +0.15 on the probability scale
        expect_equal(d$prob[d$arm == "trt2"],
                     d$pi_cluster[d$arm == "trt2"] + 0.15)
})

test_that("h1_binom logit-scale shift stays in (0, 1)", {
        set.seed(3)
        d <- h1_binom(pi = 0.3, H = 1, n = 40, rho = 0.05,
                      shift = c(0, 2), shift_scale = "logit")
        expect_true(all(d$prob > 0 & d$prob < 1))
})

test_that("h1_binom validates its inputs", {
        expect_error(h1_binom(pi = 1.2, H = 3, n = 40, rho = 0.05))
        expect_error(h1_binom(pi = 0.3, H = 3, n = 40, rho = 0))
        expect_error(h1_binom(pi = 0.3, H = 3, n = 40, rho = 1.5))
})

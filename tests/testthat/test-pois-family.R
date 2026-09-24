# h1_pois / PI_pois / SIM_pois. The EB path needs no gMAP MCMC, so these run
# fast; overdispersed historical counts keep the moment estimator off its cap.

test_that("h1_pois returns per-study/arm counts with a shared cluster rate", {
        set.seed(1)
        d <- h1_pois(lambda = 8, H = 4, offset = 3, shape = 12)
        expect_s3_class(d, "data.frame")
        expect_named(d, c("study", "arm", "events", "offset", "lambda_cluster", "rate"))
        expect_equal(nrow(d), 4)                 # control only, 4 studies
        expect_true(all(d$arm == "control"))
        expect_true(all(d$offset == 3))
})

test_that("h1_pois treatment arms share the cluster rate (rate-difference shift)", {
        set.seed(2)
        d <- h1_pois(lambda = 8, H = 1, offset = 3, shape = 12, shift = c(0, 4))
        expect_setequal(as.character(d$arm), c("control", "trt1", "trt2"))
        expect_equal(d$rate[d$arm == "control"], d$lambda_cluster[d$arm == "control"])
        expect_equal(d$rate[d$arm == "trt2"],
                     d$lambda_cluster[d$arm == "trt2"] + 4)
})

test_that("h1_pois log-scale shift gives a constant rate ratio", {
        set.seed(3)
        d <- h1_pois(lambda = 8, H = 1, offset = 3, shape = 12,
                     shift = c(0, log(2)), shift_scale = "log")
        expect_equal(unname(d$rate[d$arm == "trt2"] / d$rate[d$arm == "control"]), 2,
                     tolerance = 1e-8)
})

test_that("h1_pois validates its inputs", {
        expect_error(h1_pois(lambda = -1, H = 3, offset = 3, shape = 12))
        expect_error(h1_pois(lambda = 8, H = 3, offset = 3, shape = 0))
        expect_error(h1_pois(lambda = 8, H = 3, offset = 0, shape = 12))
})

test_that("PI_pois (EB) returns a PI object with an NA posterior", {
        histdat <- data.frame(events = c(4, 22, 6, 35, 9, 28, 3, 26), offset = rep(2, 8))
        m_pi <- suppressWarnings(PI_pois(histdat, offset0 = 3, prior_est = "EB", seed = 1))
        expect_true(inherits(m_pi, "PI") && inherits(m_pi, "MultiMAP_pois"))
        expect_named(m_pi$pred_int, c("median", "lower", "upper"))
        expect_true(is.na(m_pi$eff_n_post))
        expect_identical(m_pi$distr, "poisson")
})

test_that("SIM_pois (EB) returns a one-row data frame with the key columns", {
        out <- SIM_pois(n_sim = 8, H = 6, offset = 3, lambda = 8, shape = 12,
                        shift = c(0, 0), prior_est = "EB", contr = "rate_diff",
                        seed = 1, parallel = FALSE)
        expect_s3_class(out, "data.frame")
        expect_equal(nrow(out), 1)
        expect_true(all(c("lambda", "lambda_cur", "shape", "sd_h_between", "cv_h_between",
                          "fwer", "power", "pi_t1e", "pi_power") %in% names(out)))
        expect_false(is.na(out$fwer))
        expect_true(is.na(out$power))
})

test_that("SIM_pois converts sd_h / cv_h to the same shape", {
        lambda <- 8; shape <- 12
        o1 <- SIM_pois(n_sim = 3, H = 6, offset = 3, lambda = lambda,
                       sd_h = lambda / sqrt(shape), shift = c(0, 0),
                       prior_est = "EB", seed = 1, parallel = FALSE)
        o2 <- SIM_pois(n_sim = 3, H = 6, offset = 3, lambda = lambda,
                       cv_h = 1 / sqrt(shape), shift = c(0, 0),
                       prior_est = "EB", seed = 1, parallel = FALSE)
        expect_equal(o1$shape, shape, tolerance = 1e-8)
        expect_equal(o2$shape, shape, tolerance = 1e-8)
})

test_that("SIM_pois requires exactly one between-study specification", {
        expect_error(SIM_pois(n_sim = 2, H = 4, offset = 3, lambda = 8,
                              shift = c(0, 0), prior_est = "EB", parallel = FALSE))
        expect_error(SIM_pois(n_sim = 2, H = 4, offset = 3, lambda = 8,
                              shape = 12, sd_h = 2, shift = c(0, 0),
                              prior_est = "EB", parallel = FALSE))
})

test_that("SIM_pois conflict scenario reports pi_power, not pi_t1e", {
        out <- SIM_pois(n_sim = 8, H = 6, offset = 3, lambda = 8, shape = 12,
                        shift = c(0, 4), lambda_cur = 4, prior_est = "EB",
                        contr = "rate_diff", seed = 1, parallel = FALSE)
        expect_true(is.na(out$pi_t1e))
        expect_false(is.na(out$pi_power))
})

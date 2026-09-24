#' One-way layout data generator (Poisson endpoint)
#'
#' Simulates count data from a gamma-Poisson (negative-binomial) one-way
#' random-effects model, the Poisson analogue of \code{\link{h1_norm}} and
#' \code{\link{h1_binom}}. Each of the \code{H} studies (indexed \eqn{h = 1,
#' \ldots, H}) is a cluster with its own true rate \eqn{\lambda_h \sim
#' \mathrm{Gamma}(\mathrm{shape}, \mathrm{shape} / \lambda)}, which has mean
#' \code{lambda} and variance \eqn{\lambda^2 / \mathrm{shape}} (so the marginal
#' counts are negative binomial, following Tarone, 1982). Within a study the
#' control group (\eqn{m = 0}) has \eqn{Y_0 \sim \mathrm{Poisson}(\lambda_h \cdot
#' n_0)}; the treatment groups (\eqn{m = 1, \ldots, M}) share the same cluster
#' rate \eqn{\lambda_h} and are offset by a fixed \code{shift} (on the rate or
#' the log scale), so \code{shift_m} is the true treatment-vs-control effect on
#' that scale. The number of treatment groups \eqn{M} is defined by
#' \code{length(shift)} (\code{shift = NULL} gives controls only).
#'
#' The between-study standard deviation of the true cluster rates is
#' \eqn{\lambda / \sqrt{\mathrm{shape}}}; small \code{shape} means strong
#' heterogeneity, large \code{shape} tends to a homogeneous Poisson model.
#'
#' @param lambda grand rate (single positive value); the mean of the gamma
#'   distribution of the cluster rates.
#' @param H number of studies (clusters), each with its own drawn cluster rate.
#' @param offset exposure (number of experimental units) per group. A single
#'   value (common to all groups) or a vector giving a separate exposure per
#'   generated group (recycled to \eqn{H (M + 1)} groups), so the historical
#'   studies and the current-trial groups may have different exposures.
#' @param shape gamma shape of the between-study rate distribution (single
#'   positive value; Tarone's \code{p}). Larger \code{shape} -> less
#'   heterogeneity.
#' @param shift treatment effect(s) relative to the control, applied identically
#'   in every study on the scale given by \code{shift_scale}. A numeric vector;
#'   its length is the number of treatment groups \eqn{M} per study. \code{NULL}
#'   (default) gives controls only.
#' @param shift_scale scale on which \code{shift} acts: \code{"rate"} (default)
#'   adds it to the control rate (\code{shift} is the true rate difference),
#'   \code{"log"} adds it on the log scale (\code{exp(shift)} is the true rate
#'   ratio). Under \code{shift = 0} both scales leave the treatment arms equal to
#'   the control. On the rate scale a shift that pushes a group rate below 0 is
#'   clipped to 0, with a warning.
#'
#' @return A data frame with one row per study x group and columns \code{study}
#'   (cluster factor \code{1..H}), \code{arm} (factor: \code{"control"} or
#'   \code{"trt1"}..\code{"trt<M>"}), \code{events} (the count), \code{offset}
#'   (the exposure), \code{lambda_cluster} (the study's shared true rate) and
#'   \code{rate} (the group's true rate).
#'
#' @seealso \code{\link{SIM_pois}}, \code{\link{h1_binom}}, \code{\link{h1_norm}}
#'
#' @references
#' Tarone, R. E. (1982). The use of historical control information in testing for
#' a trend in Poisson means. \emph{Biometrics, 38}(2), 457--462.
#' \doi{10.2307/2530459}
#'
#' @examples
#' \donttest{
#' # Four historical control studies, no treatment arms
#' h1_pois(lambda = 8, H = 4, offset = 3, shape = 12)
#'
#' # One trial with a control and two treatment arms (rate-difference shift)
#' h1_pois(lambda = 8, H = 1, offset = 3, shape = 12, shift = c(0, 4))
#' }
#' @export
h1_pois <- function(lambda, H, offset, shape, shift = NULL,
                    shift_scale = c("rate", "log")){

        shift_scale <- match.arg(shift_scale)

        # Number of treatment groups M is defined by the length of shift
        M <- length(shift)

        #-----------------------------------------------------------------------
        # Input checks

        if(!is.numeric(lambda) || length(lambda) != 1 || lambda <= 0)
                stop("'lambda' must be a single positive number.")
        if(!is.numeric(shape) || length(shape) != 1 || shape <= 0)
                stop("'shape' must be a single positive number.")
        if(!is.numeric(H) || length(H) != 1 || H < 1)
                stop("'H' must be a single positive integer.")
        if(!is.numeric(offset) || any(offset <= 0))
                stop("'offset' must be positive (a single value or one per generated group).")

        #-----------------------------------------------------------------------
        # Groups per study: control (m = 0, shift 0) then the M treatment groups

        arm_labels <- c("control", if(M > 0) paste0("trt", seq_len(M)))
        arm_shift  <- c(0,         if(M > 0) shift)
        n_arm      <- 1 + M

        # One shared true rate per study (the cluster random effect):
        # Gamma(shape, rate = shape / lambda) -> mean lambda, var lambda^2 / shape.
        lambda_cluster <- rgamma(n = H, shape = shape, rate = shape / lambda)

        # Index vectors (study-major, then group)
        study_idx <- rep(x = seq_len(H),     each  = n_arm)
        arm_idx   <- rep(x = seq_len(n_arm), times = H)

        # Per-group exposures: scalar (common) or one per generated group.
        off_vec <- rep_len(offset, H * n_arm)

        lc_vec <- lambda_cluster[study_idx]
        sh_vec <- arm_shift[arm_idx]

        # Group rates: shift added on the rate scale (rate difference) or the log
        # scale (rate ratio).
        if(shift_scale == "rate"){
                raw  <- lc_vec + sh_vec
                rate <- pmax(raw, 0)
                if(any(raw < 0))
                        warning("Some shifted rates fell below 0 and were clipped to 0; ",
                                "consider shift_scale = 'log' or a smaller shift.", call. = FALSE)
        } else {
                rate <- lc_vec * exp(sh_vec)
        }

        # Number of events per group
        events <- rpois(n = length(rate), lambda = rate * off_vec)

        data.frame(study          = factor(study_idx),
                   arm            = factor(arm_labels[arm_idx], levels = arm_labels),
                   events         = events,
                   offset         = off_vec,
                   lambda_cluster = lc_vec,
                   rate           = rate)
}

#' One-way layout data generator (binomial endpoint)
#'
#' Simulates count data from a beta-binomial one-way random-effects model, the
#' binomial analogue of \code{\link{h1_norm}}. Each of the \code{H} studies
#' (indexed \eqn{h = 1, \ldots, H}) is a cluster with its own true success
#' probability \eqn{\pi_h \sim \mathrm{Beta}(a, b)}, where \eqn{a} and \eqn{b}
#' follow \code{pi} and the intra-class correlation \code{rho} through
#' \eqn{\rho = 1 / (1 + a + b)} (the \code{predint::rbbinom} parametrisation).
#' Within a study the control group (\eqn{m = 0}) has \eqn{Y_0 \sim
#' \mathrm{Bin}(n_0, \pi_h)}; the treatment groups (\eqn{m = 1, \ldots, M}) share
#' the same cluster proportion \eqn{\pi_h} and are offset by a fixed \code{shift}
#' (on the probability or logit scale), so \code{shift_m} is the true
#' treatment-vs-control effect on that scale. The number of treatment groups
#' \eqn{M} is defined by \code{length(shift)} (\code{shift = NULL} gives controls
#' only).
#'
#' The between-study variance of the true cluster proportions is
#' \eqn{\pi(1-\pi)\rho}; equivalently the between-study standard deviation is
#' \eqn{\sqrt{\pi(1-\pi)\rho}}.
#'
#' @param pi grand success probability (single value in (0, 1)); the mean of the
#'   beta distribution of the cluster proportions.
#' @param H number of studies (clusters), each with its own drawn cluster
#'   proportion.
#' @param n number of trials (patients) per group. A single value (common to all
#'   groups) or a vector giving a separate size per generated group (recycled to
#'   \eqn{H (M + 1)} groups), so the historical studies and the current-trial
#'   groups may have different sizes.
#' @param rho intra-class correlation of the beta-binomial model (single value in
#'   (0, 1)); controls the between-study variation. As \code{rho} approaches 0
#'   the model tends to a plain binomial (no heterogeneity).
#' @param shift treatment effect(s) relative to the control, applied identically
#'   in every study on the scale given by \code{shift_scale}. A numeric vector;
#'   its length is the number of treatment groups \eqn{M} per study. \code{NULL}
#'   (default) gives controls only.
#' @param shift_scale scale on which \code{shift} acts: \code{"prob"} (default)
#'   adds it to the control probability (\code{shift} is the true risk
#'   difference), \code{"logit"} adds it on the logit scale (\code{exp(shift)} is
#'   the true odds ratio). Under \code{shift = 0} both scales leave the treatment
#'   arms equal to the control. On the probability scale a shift that pushes an
#'   arm probability outside (0, 1) is clipped, with a warning.
#'
#' @return A data frame with one row per study x group and columns \code{study}
#'   (cluster factor \code{1..H}), \code{arm} (factor: \code{"control"} or
#'   \code{"trt1"}..\code{"trt<M>"}), \code{events} (number of successes),
#'   \code{n} (number of trials), \code{pi_cluster} (the study's shared true
#'   proportion) and \code{prob} (the group's true success probability).
#'
#' @seealso \code{\link{SIM_binom}}, \code{\link{h1_norm}}, \code{\link{h1_pois}}
#'
#' @examples
#' \donttest{
#' # Four historical control studies, no treatment arms
#' h1_binom(pi = 0.3, H = 4, n = 40, rho = 0.05)
#'
#' # One trial with a control and two treatment arms (risk-difference shift)
#' h1_binom(pi = 0.3, H = 1, n = 40, rho = 0.05, shift = c(0, 0.15))
#' }
#' @export
h1_binom <- function(pi, H, n, rho, shift = NULL, shift_scale = c("prob", "logit")){

        shift_scale <- match.arg(shift_scale)

        # Number of treatment groups M is defined by the length of shift
        M <- length(shift)

        #-----------------------------------------------------------------------
        # Input checks

        if(!is.numeric(pi) || length(pi) != 1 || pi <= 0 || pi >= 1)
                stop("'pi' must be a single number in (0, 1).")
        if(!is.numeric(rho) || length(rho) != 1 || rho <= 0 || rho >= 1)
                stop("'rho' must be a single number in (0, 1).")
        if(!is.numeric(H) || length(H) != 1 || H < 1)
                stop("'H' must be a single positive integer.")
        if(!is.numeric(n) || any(n < 1))
                stop("'n' must be positive (a single value or one per generated group).")

        #-----------------------------------------------------------------------
        # Beta parameters implied by (pi, rho): a + b = (1 - rho) / rho,
        # a = pi (a + b), b = (a + b) - a  (predint::rbbinom parametrisation).

        ab <- (1 - rho) / rho
        a  <- pi * ab
        b  <- ab - a

        #-----------------------------------------------------------------------
        # Groups per study: control (m = 0, shift 0) then the M treatment groups

        arm_labels <- c("control", if(M > 0) paste0("trt", seq_len(M)))
        arm_shift  <- c(0,         if(M > 0) shift)
        n_arm      <- 1 + M

        # Per-group sample sizes: scalar (common) or one per generated group.
        n_vec <- rep_len(n, H * n_arm)

        # One shared true proportion per study (the cluster random effect)
        pi_cluster <- rbeta(n = H, shape1 = a, shape2 = b)

        # Index vectors (study-major, then group)
        study_idx <- rep(x = seq_len(H),     each  = n_arm)
        arm_idx   <- rep(x = seq_len(n_arm), times = H)

        pc_vec <- pi_cluster[study_idx]
        sh_vec <- arm_shift[arm_idx]

        # Group probabilities: shift added on the probability or the logit scale.
        if(shift_scale == "prob"){
                raw  <- pc_vec + sh_vec
                prob <- pmin(pmax(raw, 0), 1)
                if(any(raw < 0 | raw > 1))
                        warning("Some shifted probabilities fell outside (0, 1) and were clipped; ",
                                "consider shift_scale = 'logit' or a smaller shift.", call. = FALSE)
        } else {
                prob <- plogis(qlogis(pc_vec) + sh_vec)
        }

        # Number of successes per group
        events <- rbinom(n = length(prob), size = n_vec, prob = prob)

        data.frame(study      = factor(study_idx),
                   arm        = factor(arm_labels[arm_idx], levels = arm_labels),
                   events     = events,
                   n          = n_vec,
                   pi_cluster = pc_vec,
                   prob       = prob)
}

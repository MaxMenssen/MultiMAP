#' One-way layout data generator (normal endpoint)
#'
#' Simulates data from a one-way random-effects model for a normally distributed
#' endpoint. For study \eqn{h = 1, \ldots, H} the control group (\eqn{m = 0}) is
#' \eqn{y_{hi} = \mu + a_h + e_{hi}} and treatment group \eqn{m} is
#' \eqn{y_{hmi} = \mu + a_h + \mathrm{shift}_m + e_{hmi}}. For every study
#' \eqn{h} a shared random effect \eqn{a_h} is drawn; the treatment groups of
#' study \eqn{h} share this same random effect and the same within-group residual
#' SD as the control, differing only by a fixed mean shift, so \code{shift_m} is
#' exactly the true treatment-vs-control effect. The number of treatment groups
#' \eqn{M} is defined by \code{length(shift)} (\code{shift = NULL} gives controls
#' only).
#'
#' @param mu grand mean (single value).
#' @param H number of studies, each defining one shared-random-effect block.
#' @param n number of observations per group (control and each treatment group).
#'   A single value (common to all groups) or a vector giving a separate size per
#'   generated group (recycled to \eqn{H (M + 1)} groups), so the historical
#'   studies and the current-trial groups may have different sizes.
#' @param sd_h between-study (random effect) standard deviation (single value).
#' @param sd_hi within-group residual standard deviation. Either a single value
#'   (common residual SD for all studies) or a vector of length \code{H} giving a
#'   separate residual SD per study. Each study's treatment groups inherit its
#'   residual SD.
#' @param shift mean shift(s) of the treatment groups relative to their control. A
#'   numeric vector; its length is the number of treatment groups \eqn{M} per
#'   study. \code{NULL} (default) gives no treatment groups. Applied identically
#'   across all studies.
#'
#' @return A data frame with columns \code{y_hi} (the observation), \code{a}
#'   (study factor \code{1..H}, the shared random effect / cluster), \code{arm}
#'   (factor: \code{"control"} or \code{"trt1"}..\code{"trt<M>"}) and \code{e}
#'   (unique observation id).
#'
#' @seealso \code{\link{SIM_norm}}, \code{\link{h1_binom}}, \code{\link{h1_pois}}
#'
#' @examples
#' \donttest{
#' dat <- h1_norm(mu = 10, H = 5, n = 30, sd_h = 2, sd_hi = 3,
#'                shift = c(0, 1.5))
#' head(dat)
#' }
#' @export
h1_norm <- function(mu, H, n, sd_h, sd_hi, shift = NULL){

        # Number of treatment groups M is defined by the length of shift
        M <- length(shift)

        # Residual SD per study: recycle a single value to one per study,
        # otherwise require exactly one value per study.
        if(length(sd_hi) == 1){
                sd_hi <- rep(x = sd_hi, times = H)
        }
        if(length(sd_hi) != H){
                stop("sd_hi must have length 1 or H (one residual SD per study).")
        }

        # Groups per study: control (m = 0, shift 0) then the M treatment groups
        arm_labels <- c("control", if(M > 0) paste0("trt", seq_len(M)))
        arm_shift  <- c(0,         if(M > 0) shift)
        n_arm      <- 1 + M

        # Per-group sample sizes: scalar (common) or one per generated group.
        n_groups <- H * n_arm
        n_vec    <- rep_len(n, n_groups)

        # One shared random effect per study
        a_h <- rnorm(n = H, sd = sd_h)

        # Group table (study-major, then group), then expand to observations.
        grp_study <- rep(x = seq_len(H),     each  = n_arm)
        grp_arm   <- rep(x = seq_len(n_arm), times = H)
        obs_group <- rep(seq_len(n_groups), times = n_vec)

        # General mean + shared random effect + treatment shift; residual SD by study
        mean_vec <- mu + a_h[grp_study[obs_group]] + arm_shift[grp_arm[obs_group]]
        e_hi     <- rnorm(n = length(obs_group), sd = sd_hi[grp_study[obs_group]])
        y_hi     <- mean_vec + e_hi

        data.frame(y_hi = y_hi,
                   a    = factor(grp_study[obs_group]),
                   arm  = factor(arm_labels[grp_arm[obs_group]], levels = arm_labels),
                   e    = factor(seq_along(y_hi)))
}

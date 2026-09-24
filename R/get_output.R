#' Extract a component from a MultiMAP object
#'
#' Single entry point that returns a chosen component of a MultiMAP object -- the
#' output of \code{\link{SCI_norm}}, \code{\link{PI_norm}}, \code{\link{SCI_binom}},
#' \code{\link{PI_binom}} or \code{\link{SCI_pois}} -- selected via \code{which}.
#' Not every component is populated by every constructor: the \code{PI_*}
#' functions leave the contrast-related components (\code{"SCI"},
#' \code{"contr_mat"}, \code{"post"}) as \code{NA}.
#'
#' @param x a MultiMAP object (from \code{\link{SCI_norm}}, \code{\link{PI_norm}},
#'   \code{\link{SCI_binom}}, \code{\link{PI_binom}} or \code{\link{SCI_pois}}).
#' @param which which component to extract:
#'   \describe{
#'     \item{\code{"SCI"}}{simultaneous credible intervals of the contrasts
#'       (\code{x$comp}); the default.}
#'     \item{\code{"weights"}}{prior and posterior mixture weights (\code{x$weights}).}
#'     \item{\code{"ess"}}{effective sample sizes of the (non-robustified) prior,
#'       the robustified prior and the concurrent-control posterior, named
#'       \code{"MAP ESS"} / \code{"Robust. MAP ESS"} / \code{"Post. ESS"}.}
#'     \item{\code{"contr_mat"}}{the contrast matrix (\code{x$contr_mat}).}
#'     \item{\code{"prior"}}{list of the prior and its robustified version,
#'       named \code{"Prior"} / \code{"Prior (rob.)"}.}
#'     \item{\code{"post"}}{the concurrent-control posterior mixture
#'       (\code{x$posterior}).}
#'     \item{\code{"prior_pred"}}{the (prior) predictive distribution
#'       (\code{x$pred_dist}).}
#'     \item{\code{"pred_int"}}{the prediction interval (\code{x$pred_int}).}
#'   }
#'
#' @return The requested component.
#'
#' @examples
#' \donttest{
#' histdat <- data.frame(mean = c(9.8, 10.2, 9.5, 10.0),
#'                       sd   = c(2.1, 1.9, 2.3, 2.0),
#'                       n    = c(40, 55, 38, 47))
#' newdat  <- data.frame(group = c("control", "d1", "d454"),
#'                       mean  = c(10.1, 11.4, 12.7),
#'                       sd    = c(2.0, 2.2, 2.1),
#'                       n     = c(30, 30, 30))
#' res <- SCI_norm(histdat, newdat, contr = "mean_diff", seed = 1)
#' get_output(res)                     # SCI (default)
#' get_output(res, which = "weights")
#' get_output(res, which = "pred_int")
#' }
#' @export
get_output <- function(x,
                         which = c("SCI", "weights", "ess", "contr_mat",
                                   "prior", "post", "prior_pred", "pred_int")){

        # Stop if x is not a MultiMAP object
        if(!inherits(x, "MultiMAP")){
                stop("x is not of class MultiMAP")
        }

        which <- match.arg(which)

        # Multi-arm objects (histdat supplied as a named list) carry the borrowing
        # components per arm; the accessors then return per-arm structures.
        multi <- !is.null(x$borrow_arms)

        if(multi){
                arms <- x$borrow_arms
                return(switch(which,
                       SCI        = x$comp,
                       weights    = x$weights,           # named list of data frames
                       ess        = {
                               m <- cbind("MAP ESS"         = x$eff_n_map,
                                          "Robust. MAP ESS" = x$eff_n,
                                          "Post. ESS"       = x$eff_n_post)
                               rownames(m) <- arms
                               as.data.frame(m)
                       },
                       contr_mat  = x$contr_mat,
                       prior      = stats::setNames(
                                       lapply(arms, function(a)
                                              stats::setNames(list(x$prior[[a]], x$prior_rob[[a]]),
                                                              c("Prior", "Prior (rob.)"))),
                                       arms),
                       post       = x$posterior,         # named list of posteriors
                       prior_pred = x$pred_dist,         # named list
                       pred_int   = x$pred_int))         # named list
        }

        switch(which,
               SCI        = x$comp,
               weights    = x$weights,
               ess        = {
                       ess <- c(x$eff_n_map, x$eff_n, x$eff_n_post)
                       names(ess) <- c("MAP ESS", "Robust. MAP ESS", "Post. ESS")
                       ess
               },
               contr_mat  = x$contr_mat,
               prior      = {
                       p_list <- list(x$prior, x$prior_rob)
                       names(p_list) <- c("Prior", "Prior (rob.)")
                       p_list
               },
               post       = x$posterior,
               prior_pred = x$pred_dist,
               pred_int   = x$pred_int)
}

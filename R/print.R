#' Print a MultiMAP object
#'
#' Compact console output for a MultiMAP object (normal or binomial endpoint).
#' For \code{SCI} output (class \code{"SCI"}, from \code{\link{SCI_norm}} /
#' \code{\link{SCI_binom}}) it prints the table of simultaneous credible
#' intervals (\code{x$comp}); for \code{PI} output (class \code{"PI"}, from
#' \code{\link{PI_norm}} / \code{\link{PI_binom}}) it prints the prediction
#' interval (\code{x$pred_int}).
#'
#' @param x a MultiMAP object (from \code{\link{SCI_norm}},
#'   \code{\link{PI_norm}}, \code{\link{SCI_binom}}, \code{\link{PI_binom}} or
#'   \code{\link{SCI_pois}}).
#' @param ... passed to the underlying \code{print} call.
#'
#' @return The object \code{x}, invisibly.
#'
#' @examples
#' \donttest{
#' histdat <- data.frame(mean = c(9.8, 10.2, 9.5, 10.0),
#'                       sd   = c(2.1, 1.9, 2.3, 2.0),
#'                       n    = c(40, 55, 38, 47))
#' pred <- PI_norm(histdat = histdat, n0 = 30, seed = 1)
#' print(pred)
#' }
#' @export
print.MultiMAP <- function(x, ...){

        # Stop if x is not a MultiMAP object
        if(!inherits(x, "MultiMAP")){
                stop("x is not of class MultiMAP")
        }

        if(inherits(x, "SCI")){
                # SCI_*(): simultaneous credible intervals (one table, any # arms)
                print(x$comp, ...)
        } else if(!is.null(x$borrow_arms)){
                # Multi-arm PI_*(): one prediction-interval row per borrowed arm
                print(as.data.frame(do.call(rbind, x$pred_int)), ...)
        } else {
                # PI_*(): single prediction interval
                print(x$pred_int, ...)
        }

        invisible(x)
}

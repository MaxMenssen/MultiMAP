#' Summarise a MultiMAP object
#'
#' Printed summary of a MultiMAP object (normal, binomial or Poisson endpoint).
#' The borrowing information is printed as one block per arm for which historical
#' data were used -- the concurrent control alone when \code{histdat} is a single
#' data frame, or every named arm when \code{histdat} is a list (multi-arm
#' borrowing). Each block reports, at the \eqn{100(1-\alpha)\%} level, the
#' prediction interval (\code{object$pred_int}) for that arm's mean (normal),
#' proportion (binomial) or rate (Poisson), and its effective sample sizes (the
#' MAP prior, the robustified prior and -- for \code{\link{SCI_norm}} /
#' \code{\link{SCI_binom}} / \code{\link{SCI_pois}} output, class \code{"SCI"} --
#' the posterior). For \code{"SCI"} output the observed value of the arm and
#' whether it is covered by the prediction interval are added, and the
#' simultaneous credible intervals of the contrasts (\code{object$comp}) are
#' printed after the blocks, with the heading naming the contrast type taken from
#' \code{object$contr}.
#'
#' @param object a MultiMAP object (from \code{\link{SCI_norm}},
#'   \code{\link{PI_norm}}, \code{\link{SCI_binom}}, \code{\link{PI_binom}} or
#'   \code{\link{SCI_pois}}).
#' @param ... currently ignored.
#'
#' @return The object, invisibly.
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
#' summary(res)
#' }
#' @export
summary.MultiMAP <- function(object, ...){

        # Stop if object is not a MultiMAP object
        if(!inherits(object, "MultiMAP")){
                stop("object is not of class MultiMAP")
        }

        # Coverage level 1 - alpha (falls back to the literal "1-alpha" if the
        # object carries no alpha).
        conf <- if(is.null(object$alpha)) "1-alpha" else 1 - object$alpha

        # The concurrent control location is a mean (normal), a proportion
        # (binomial) or a rate (Poisson); name it accordingly.
        unit <- switch(as.character(object$distr),
                       betabinomial = "proportion",
                       poisson      = "rate",
                       "mean")

        # SCI_*() objects carry class "SCI" (observed value + coverage + posterior
        # ESS); PI_*() objects do not.
        is_sci <- inherits(object, "SCI")

        # One borrowing block per arm, identical layout for single- and multi-arm
        # output: the prediction interval, the observed value + coverage (SCI
        # only) and the effective sample sizes.
        arm_block <- function(label, pit, obs, e_map, e_rob, e_post){
                cat("\nArm '", label, "'\n", sep = "")
                cat("  ", 100*conf, " % prediction interval for the ", unit, ": [",
                    signif(pit["lower"], 4), ", ", signif(pit["upper"], 4), "]\n", sep = "")
                if(is_sci){
                        lwr <- if(is.na(pit["lower"])) -Inf else pit["lower"]
                        upr <- if(is.na(pit["upper"]))  Inf else pit["upper"]
                        covered <- obs >= lwr & obs <= upr
                        cat("  observed ", unit, ": ", signif(obs, 4),
                            if(covered) " (covered)" else " (NOT covered)", "\n", sep = "")
                }
                cat("  ESS  MAP: ", e_map, " | robust. MAP: ", e_rob,
                    if(is_sci) paste0(" | posterior: ", e_post) else "", "\n", sep = "")
        }

        # Simultaneous credible intervals of the contrasts (SCI output only).
        print_sci <- function(){
                if(is_sci){
                        cat("\n", 100*conf, " % simultaneous credible intervals for ",
                            .contr_label(object$contr), "\n", sep = "")
                        print(object$comp)
                }
        }

        #---- multi-arm objects: one borrowing block per borrowed arm ----------
        if(!is.null(object$borrow_arms)){
                arms <- object$borrow_arms
                cat("Multi-arm borrowing for ", length(arms), " arm(s): ",
                    paste(arms, collapse = ", "), "\n", sep = "")
                for(a in arms){
                        obs <- if(is_sci)
                                object$newdat$mean[match(a, as.character(object$newdat$group))]
                               else NA
                        arm_block(a, object$pred_int[[a]], obs,
                                  object$eff_n_map[a], object$eff_n[a], object$eff_n_post[a])
                }
                print_sci()
                return(invisible(object))
        }

        #---- single-arm objects: one block for the concurrent control ---------
        # The borrowed arm is the concurrent control (row 1 of newdat for SCI
        # output; PI output carries no newdat).
        cat("Single-arm borrowing for the concurrent control\n", sep = "")
        label <- if(is_sci) as.character(object$newdat$group[1]) else "concurrent control"
        obs   <- if(is_sci) object$newdat$mean[1] else NA
        arm_block(label, object$pred_int, obs,
                  object$eff_n_map, object$eff_n, object$eff_n_post)
        print_sci()

        invisible(object)
}

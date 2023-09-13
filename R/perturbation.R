# ------------------------------------------------------------------------------
# Perturbation bootstrap for standard errors
# ------------------------------------------------------------------------------

#' Perturbation bootstrap for calculating the standard errors of the ROC estimates
#'
#' @param nbt number of replicates for perturbation boostrap
#' @param S phenotyping score S for all, including labeled and unlabeled set
#' @param Y outcome Y for all, including labeled and unlabeled set; containing Y NA
#' @param method validation methods, choosing from supROC and ssROC
#' @return list containing ROC estimates for each perturbation bootstrap
#' @export

perturbation <- function(nbt, S, Y, method) {
  res_bt <- NULL

  # get labeled data
  labeled_ind <- which(!is.na(Y))
  Y_labeled <- Y[labeled_ind]
  S_labeled <- S[labeled_ind]
  
  for (ibt in 1:nbt) {
    ptb_wgt <- 4 * rbeta(length(S_labeled), 1 / 2, 3 / 2)
    if (method == "supROC") {
      res_bt[[ibt]] <- tryCatch(
        supROC(S_labeled, Y_labeled, W_labeled = ptb_wgt),
        error = function(e) NA
      )
    } else if (method == "ssROC") {
      res_bt[[ibt]] <- tryCatch(
        ssROC(S, Y, W_labeled = ptb_wgt),
        error = function(e) NA
      )
    }
  }
  return(res_bt)
}

# -----------------------------------------------------------------------------
# Relative Efficiency
# -----------------------------------------------------------------------------

#' Compute relative efficiency of supROC vs ssROC
#' @param sl_pert perturbation result of supROC
#' @param ssl_pert pertubation result of ssROC
#' @param fpr_threshold FPR threshold to compute ROC parameters. Must choose 
#' from seq(0.01, 0.99, by = 0.01)
#' @return vector of relative efficiency computed at each ROC parameter
#' @export

re <- function(sl_pert, ssl_pert, fpr_threshold) {
  if (!(fpr_threshold %in% round(seq(0.01, 0.99, by = 0.01), 2))) {
    stop("Input fpr_threshold is not allowed.")
  }
  if (length(sl_pert) < 500) {
    warning("Number of boostrap too small, RE estimates are not accurate.
            Increase the number of bootstrap")
  }
  sl <- do.call(rbind, lapply(
    sl_pert,
    function(ll) ll[round(ll[, "FPR"], 2) == fpr_threshold, ]
  ))
  ssl <- do.call(rbind, lapply(
    ssl_pert,
    function(ll) ll[round(ll[, "FPR"], 2) == fpr_threshold, ]
  ))
  out <- apply(sl, 2, var) / apply(ssl, 2, var)
  return(out[!names(out) == "FPR"])
}

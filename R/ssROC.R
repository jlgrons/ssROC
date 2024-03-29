# -----------------------------------------------------------------------------
# SS-ROC Imputation Method
# -----------------------------------------------------------------------------
#'
#' Imputation based semi-supervised method
#' @param S phenotyping score S for all, including labeled and unlabeled set
#' @param Y outcome Y for all, including labeled and unlabeled set; containing Y NA
#' @param W_labeled optional vector of weights for labeled set
#' @param W_unlabeled optional vector of weights for unlabeled set
#' @param fpr_vals desired fpr sequence for output
#' @param bandwidth bandwidth for smoothing
#' @param transform_score Whether to use ecdf to transform score
#' @return matrix containing ssROC estimates at fpr_vals
#' @export

ssROC <- function(S, Y,
                  W_labeled = NULL,
                  W_unlabeled = NULL,
                  fpr_vals = seq(0.01, 0.99, by = 0.01),
                  bandwidth = NULL,
                  transform_score = TRUE) {
  id_labeled <- which(!is.na(Y))

  Y_labeled <- Y[id_labeled]

  if(transform_score){
    ecdf_S <- stats::ecdf(S)
    S <- ecdf_S(S)
  }

  S_labeled <- S[id_labeled]
  S_unlabeled <- S[-id_labeled]

  N_unlabeled <- length(S_unlabeled)
  n_labeled <- length(Y_labeled)

  if (is.null(W_labeled)) {
    W_labeled <- rep(1, n_labeled)
  }

  if (is.null(W_unlabeled)) {
    W_unlabeled <- rep(1, N_unlabeled)
  }

  if (is.null(bandwidth)) {
    bandwidth <- sqrt(Hmisc::wtd.var(S_labeled, W_labeled)) / (n_labeled^0.45)
  }

  mhat <- npreg(S_labeled, Y_labeled, S_unlabeled, bandwidth, Wt = W_labeled)

  result <- roc(
    S = S_unlabeled, Y = mhat, W = W_unlabeled, fpr_vals = fpr_vals
  )

  return(result)
}

# -----------------------------------------------------------------------------
# Supervised validation
# -----------------------------------------------------------------------------

#' ROC analysis using supervised method
#'
#' @param S_labeled phenotyping score S in labeled set
#' @param Y_labeled outcome Y in labeled set
#' @param W_labeled optional vector of weights for labeled set
#' @param fpr_vals desired fpr sequence for output
#' @param ecdf_transform Whether conducting cdf transformation
#' @return List containing AUC and ROC estimates
#' @importFrom stats ecdf
#' @export

supervised <- function(S_labeled, Y_labeled, W_labeled = NULL,
                           fpr_vals = seq(0.01, 0.99, by = 0.01),
                           ecdf_transform = TRUE){

  if(is.null(W_labeled)){

    W_labeled <- rep(1, length(S_labeled))

    }

  if(ecdf_transform){

    ecdf_S <- ecdf(S_labeled)
    S_labeled <- ecdf_S(S_labeled)

  }


  result <- interpolated_ROC(S = S_labeled, Y = Y_labeled, W = W_labeled,
                   fpr_vals = fpr_vals)

  return(result)

}

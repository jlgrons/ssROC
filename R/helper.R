# -----------------------------------------------------------------------------
# Helper Functions
# -----------------------------------------------------------------------------

#' Define logit function
#' @param x a numeric object
#' @noRd

logit <- function(x) {
  log(x / (1 - x))
}

#' Define expit function
#' @param x a numeric object
#' @noRd

expit <- function(x) {
  1 / (1 + exp(-x))
}

#' Nonparametric local constant regression
#' @param S_t phenotyping score S in the labeled dataset
#' @param Y_t outcome Y in the labeled dataset
#' @param S_v phenotyping score S in the unlabeled dataset
#' @param bw bandwidth for kernel smoothing
#' @param Wt optional weights
#' @param kern.mat optional kernel matrix
#' @return vector of predicted value of S_v
#' @noRd

npreg <- function(St, Yt, Sv, bw, Wt = NULL, kern.mat = NULL) {
  nv <- length(Sv)
  nt <- length(St)
  if (is.null(Wt)) {
    Wt <- rep(1, nt)
  }
  if (is.null(kern.mat)) {
    kern.mat <- sapply(1:nt, function(kk) dnorm(Sv - rep(St[kk], nv), sd = bw))
  }
  nw.est <- kern.mat %*% (Yt * Wt) * (1 / (kern.mat %*% Wt))
  return(nw.est)
}

#' Compute sum of indicator function
#' return sum_c I(Y_i FUN c)*V_i
#' @param yy Vector of cuts c
#' @param FUN relation between yy and Yi
#' @param Yi Observed outcome
#' @param Vi target vector
#' @return Vector of sum computed at each yy
#' @noRd

sum.I <- function(yy, FUN, Yi, Vi = NULL) {
  if (FUN == "<" | FUN == ">=") {
    yy <- -yy
    Yi <- -Yi
  }
  # for each distinct ordered failure time t[j], number of Xi < t[j]
  pos <- rank(c(yy, Yi), ties.method = "f")[1:length(yy)] -
    rank(yy, ties.method = "f")
  if (substring(FUN, 2, 2) == "=") pos <- length(Yi) - pos # number of Xi>= t[j]
  if (!is.null(Vi)) {
    # if FUN contains '=', tmpind is the order of descending
    if (substring(FUN, 2, 2) == "=") tmpind <- order(-Yi) else tmpind <- order(Yi)
    Vi <- apply(as.matrix(Vi)[tmpind, , drop = F], 2, cumsum)
    return(rbind(0, Vi)[pos + 1, ])
  } else {
    return(pos)
  }
}

#' Obtain ROC estimates
#' @param S Phenotyping score S
#' @param Y Outcome variable, can be imputed value
#' @param W Optional weight
#' @param fpr_vals Optional FPR levels
#' @return Matrix of ROC paremeters at each fpr_vals
#' @noRd

roc <- function(S, Y, W = NULL, fpr_vals = seq(0.01, 0.99, by = 0.01)) {
  cuts <- unique(sort(S))
  nc <- length(cuts)

  TPR_c <- sum.I(cuts, "<=", S, Y * W) / sum(Y * W)
  FPR_c <- sum.I(cuts, "<=", S, (1 - Y) * W) / sum((1 - Y) * W)

  auc <- sum(TPR_c[-1] * (FPR_c[-nc] - FPR_c[-1]))

  mu1 <- sum(Y * W) / sum(W)
  mu0 <- 1 - mu1

  PPV_c <- (TPR_c * mu1) / (TPR_c * mu1 + FPR_c * mu0)
  NPV_c <- ((1 - FPR_c) * mu0) / ((1 - FPR_c) * mu0 + (1 - TPR_c) * mu1)

  roc_c <- cbind(
    "Threshold" = cuts, "FPR" = FPR_c, "TPR" = TPR_c,
    "PPV" = PPV_c, "NPV" = NPV_c
  )

  roc <- sapply(1:ncol(roc_c), function(kk) {
    approx(roc_c[, "FPR"], roc_c[, kk],
      fpr_vals,
      rule = 2
    )$y
  })

  roc_all <- cbind(roc, auc)
  colnames(roc_all) <- c(colnames(roc_c), "AUC")
  return(roc_all)
}
# -----------------------------------------------------------------------------
  # Data generation
  # -----------------------------------------------------------------------------

#' Generating simulated data sets
#'
#' @param n_labeled the size of the labeled set
#' @param N_labeled the size of the unlabeled set
#' @param setting the simulation setting
#' @return a simulated datasets containing outcome Y and S
#' @importFrom stats ecdf
#' @export

data_generation <- function(n_labeled, N_unlabeled,setting = 1){

  N_total <- n_labeled + N_unlabeled

  if(setting == 1){

    # low AUC
    p <- 10
    rho <- 0.2
    Sigma0 <- 3*(rho+(1-rho)*diag(p));
    b0 <- c(-2, 0.1, 0.1, 0.2, -0.2, rep(0, p - 4))
    
    X <- MASS::mvrnorm(N_total, rep(0,10), Sigma0) + rbinom(N_total, 3, 0.3)
    Y <- rbinom(N_total, 1, expit(cbind(1, X, X[,3]*X[,4]) %*% c(b0, 0.1)))
    S <- expit(cbind(1, X) %*% b0)
    
    ecdf_S <- ecdf(S)
    S <- ecdf_S(S)

  }else if(setting == 2){

    # high auc
    p <- 10
    rho <- 0.2
    Sigma0 <- 3*(rho+(1-rho)*diag(p));
    b0 <- c(-4, 1, 1, 0.5, 0.5, rep(0, p - 4))
    X <- MASS::mvrnorm(N_total, rep(0,10), Sigma0) + rbinom(N_total, 3, 0.3)
    Y <- rbinom(N_total, 1, expit(cbind(1, X) %*% b0))
    S <- expit(cbind(1, X) %*% b0)

    ecdf_S <- ecdf(S)
    S <- ecdf_S(S)

  }

  Y_miss <- Y
  Y_miss[sample(c(1:N_total), N_unlabeled, replace = F)] <- NA

  my_data <- cbind(Y = Y, S = S, Y_miss = Y_miss)
  my_data <- data.frame(my_data)

  return(my_data)
}


expit <- function(x){1/(1+exp(-x))}

library(MASS)

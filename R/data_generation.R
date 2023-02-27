# -----------------------------------------------------------------------------
  # Data generation
  # -----------------------------------------------------------------------------

#' Generating simulated data sets
#'
#' @param n_labeled the size of the labeled set
#' @param N_labeled the size of the unlabeled set
#' @param mean_1 mean for controls
#' @param sd_1 sd for controls
#' @param mean_2 mean for cases
#' @param sd_2 sd for cases
#' @param prevalence prevalence of cases
#' @param misspeci whether the model is correctly specified
#' @return a simulated datasets containing outcome Y and phenotyping score for both labeled and unlabele sets
#' @importFrom stats ecdf
#' @export

data_generation <- function(n_labeled, N_unlabeled, mean_1, sd_1,
                            mean_2, sd_2, prevalence, setting = 1){

  N_total <- n_labeled + N_unlabeled

  # Y <- rbinom(N_total, size = 1, p = prevalence)

  # prop <- 0.5

  if(setting == 1){

    # low AUC
    p <- 10
    rho <- 0.2
    Sigma0 <- 3*(rho+(1-rho)*diag(p));
    b0 <- c(-2, 0, 0, 0.2, 0.2, rep(0, p - 4))

    X <- mvrnorm(N_total, rep(0,10), Sigma0) + rbinom(N_total, 3, 0.3)
    Y <- rbinom(N_total, 1, expit(cbind(1, X, X[,3]*X[,4]) %*% c(b0, 1)))
    S <- expit(cbind(1, X) %*% b0)

    # low auc
    # S <- Y * (expit(c(rep(1, floor(N_total * prop)), rep(0, ceiling(N_total * prop))) * rnorm(N_total, 4, 3)  +
    #                   c(rep(0, floor(N_total * prop)), rep(1, ceiling(N_total * prop))) * rnorm(N_total, 2, 1))) +
    #   (1 - Y) * (expit(c(rep(1, floor(N_total * prop)), rep(0, ceiling(N_total * prop))) * rnorm(N_total, 0, 1)  +
    #                      c(rep(0, floor(N_total * prop)), rep(1, ceiling(N_total * prop))) * rnorm(N_total, 2, 2)))

    # plot(density(S))
    # plot(density(S[Y == 1]), xlim = c(-0, 1))
    # lines(density(S[Y == 0]), col = "red")
    # roc_oracle <- supervised(S, Y)
    # roc_oracle[1, "AUC"]
    # ReliabilityDiagram(S, Y, plot = T)

    ecdf_S <- ecdf(S)
    S <- ecdf_S(S)

  }else if(setting == 2){

    # high auc
    p <- 10
    rho <- 0.2
    Sigma0 <- 3*(rho+(1-rho)*diag(p));
    b0 <- c(-4, 1, 1, 0.5, 0.5, rep(0, p - 4))
    X <- mvrnorm(N_total, rep(0,10), Sigma0) + rbinom(N_total, 3, 0.3)
    Y <- rbinom(N_total, 1, expit(cbind(1, X) %*% b0))
    S <- expit(cbind(1, X) %*% b0)

    # S <- Y * (expit(c(rep(1, floor(N_total * prop)), rep(0, ceiling(N_total * prop))) * rnorm(N_total, 4, 0.5)  +
    #                   c(rep(0, floor(N_total * prop)), rep(1, ceiling(N_total * prop))) * rnorm(N_total, 2, 1))) +
    #   (1 - Y) * (expit(c(rep(1, floor(N_total * prop)), rep(0, ceiling(N_total * prop))) * rnorm(N_total, -0.5, 0.5)  +
    #                      c(rep(0, floor(N_total * prop)), rep(1, ceiling(N_total * prop))) * rnorm(N_total, 1, 1)))

    # plot(density(S))
    # plot(density(S[Y == 1]), xlim = c(-0, 1))
    # lines(density(S[Y == 0]), col = "red")
    # roc_oracle <- supervised(S, Y)
    # roc_oracle[1, "AUC"]
    # ReliabilityDiagram(S, Y, plot = T)

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

# -----------------------------------------------------------------------------
# Data generation
# -----------------------------------------------------------------------------

#' Generating simulated data sets
#'
#' @param n_labeled the size of the labeled set
#' @param N_labeled the size of the unlabeled set
#' @param setting the simulation setting
#' @return a simulated datasets containing outcome Y and S
#' @export

data_generation <- function(n_labeled, N_unlabeled, setting) {
  N_total <- n_labeled + N_unlabeled

  Y <- rbinom(N_total, 1, 0.3)

  if (setting == "calibrated_highauc") {
    # perfect calibration; high AUC
    mu_1 <- 0.5
    mu_0 <- -0.5
    sigma <- 0.5
    S_uncalibrated <- Y * (rnorm(N_total, mu_1, sigma)) +
      (1 - Y) * (rnorm(N_total, mu_0, sigma))
    # Calibrated
    gamma <- (mu_1 - mu_0)/(sigma)^2
    m <- (mu_1 + mu_0)/2
    S <- 1/(1 + exp(-gamma * (S_uncalibrated - m)) * (0.7/0.3))
  } else if (setting == "calibrated_lowauc") {
    # perfect calibration; low AUC
    mu_1 <- 0.25
    mu_0 <- -0.25
    sigma <- 0.5
    S_uncalibrated <- Y * (rnorm(N_total, mu_1, sigma)) +
      (1 - Y) * (rnorm(N_total, mu_0, sigma))
    # Calibrated
    gamma <- (mu_1 - mu_0)/(sigma)^2
    m <- (mu_1 + mu_0)/2
    S <- 1/(1 + exp(-gamma * (S_uncalibrated - m)) * (0.7/0.3))
  } else if (setting == "over_highauc"){
    # over-estimation of risk; high AUC
    S <- expit(Y * (rnorm(N_total, 2.3, 0.5) + rbinom(N_total, 1, 0.3)) +
                 (1 - Y) * (rnorm(N_total, 1, 0.5) + rbinom(N_total, 1, 0.3)))
  } else if (setting == "under_highauc"){
    # under-estimation of risk; low AUC
    S <- expit(Y * (rnorm(N_total, -1.5, 0.5) + rbinom(N_total, 1, 0.1)) +
                 (1 - Y) * (rnorm(N_total, -2.6, 0.5) + rbinom(N_total, 1, 0.1)))
  } else if (setting == "over_lowauc"){
    # over-estimation of risk; low AUC
    S <- expit(Y * (rnorm(N_total, 1.2, 0.5) + rbinom(N_total, 1, 0.5)) +
                 (1 - Y) * (rnorm(N_total, 0.5, 0.5) + rbinom(N_total, 1, 0.5)))
  } else if (setting == "under_lowauc"){
    # under-estimation of risk; low AUC
    S <- expit(Y * (rnorm(N_total, -1.5, 1) + rbinom(N_total, 1, 0.3)) +
                 (1 - Y) * (rnorm(N_total, -2.5, 1) + rbinom(N_total, 1, 0.3)))
  } else if (setting == "independent"){
    # Generate S from setting 1
    mu_1 <- 0.5
    mu_0 <- -0.5
    sigma <- 0.5
    S_uncalibrated <- Y * (rnorm(N_total, mu_1, sigma)) +
      (1 - Y) * (rnorm(N_total, mu_0, sigma))
    gamma <- (mu_1 - mu_0)/(sigma)^2
    m <- (mu_1 + mu_0)/2
    S <- 1/(1 + exp(-gamma * (S_uncalibrated - m)) * (0.7/0.3))
    # random sample
    S <- sample(S)
  }

  Y_miss <- Y
  Y_miss[sample(c(1:N_total), N_unlabeled, replace = F)] <- NA

  my_data <- cbind(Y = Y, S = S, Y_miss = Y_miss)
  my_data <- data.frame(my_data)

  return(my_data)
}


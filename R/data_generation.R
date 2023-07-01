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

data_generation <- function(n_labeled, N_unlabeled, setting = 1) {
  N_total <- n_labeled + N_unlabeled

  Y <- rbinom(N_total, 1, 0.3)

  if (setting == 1) {
    # low AUC
    S <- expit(
      Y * (rnorm(N_total, 0, 0.5) + rbinom(N_total, 2, 0.3)) +
        ((1 - Y) * rnorm(N_total, 0, 0.5) + rbinom(N_total, 1, 0.3))
    )
  } else if (setting == 2) {
    # high AUC
    S <- expit(
      Y * (rnorm(N_total, 0.15, .25) + rbinom(N_total, 1, 0.3)) +
        (1 - Y) * (rnorm(N_total, 0, 0.25) + rbinom(N_total, 1, 0.1))
    )
  }

  ecdf_S <- stats::ecdf(S)
  S <- ecdf_S(S)

  Y_miss <- Y
  Y_miss[sample(c(1:N_total), N_unlabeled, replace = F)] <- NA

  my_data <- cbind(Y = Y, S = S, Y_miss = Y_miss)
  my_data <- data.frame(my_data)

  return(my_data)
}
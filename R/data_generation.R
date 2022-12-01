-----------------------------------------------------------------------------
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
#' @export

data_generation <- function(n_labeled, N_unlabeled, mean_1, sd_1,
                            mean_2, sd_2, prevalence,
                            misspec = FALSE){

  N_total <- n_labeled + N_unlabeled

  Y <- rbinom(N_total, size=1, p=prevalence)

  if(!misspec){
    S = Y*rnorm(N_total, mean_1, sd_1) +
                (1-Y)*rnorm(N_total, mean_2, sd_2)
  }else if(misspec == 1){
    S = Y * (expit(c(rep(1, floor(N_total/2)), rep(0, ceiling(N_total/2))) * rnorm(N_total, 6, 2)  +
                         c(rep(0, floor(N_total/2)), rep(1, ceiling(N_total/2))) * rnorm(N_total, 2, 0.5))) +
      (1 - Y) * (expit(c(rep(1, floor(N_total/2)), rep(0, ceiling(N_total/2))) * rnorm(N_total, 2, 2)  +
                                 c(rep(0, floor(N_total/2)), rep(1, ceiling(N_total/2))) * rnorm(N_total, 1, 0.5)))
  }else if(misspec == 2){

    S = Y*  (clog(c(rep(1, floor(N_total/2)), rep(0, ceiling(N_total/2))) * rnorm(N_total, 1.5, 2)  +

            c(rep(0, floor(N_total/2)), rep(1, ceiling(N_total/2))) * rnorm(N_total, 1, 1))) +

      (1 - Y) * (expit(c(rep(1, floor(N_total/2)), rep(0, ceiling(N_total/2))) * rnorm(N_total, 3, 2)  +

                                 c(rep(0, floor(N_total/2)), rep(1, ceiling(N_total/2))) * rnorm(N_total, -1, 0.5)))
  }

  Y_miss = Y
  Y_miss[sample(c(1:N_total), N_unlabeled, replace = F)] = NA

  my_data <- cbind(Y = Y, S = S, Y_miss = Y_miss)
  my_data = data.frame(my_data)

  return(my_data)
}


clog <- function(x){1 - exp(-exp(x))}

# ssROC

R package for Semi-Supervised ROC (ssROC) Analysis for Reliable and Streamlined Evaluation of Phenotyping Algorithms.


# Installation

```{R, eval = FALSE}
devtools::install_github(repo = "https://github.com/jlgrons/ssROC")
```

# Example
```{R, eval = FALSE}
library(ssROC)

## Set up parameters.
n_labeled <- 200
N_unlabeled <- 10000
setting <- "calibrated_highauc"

set.seed(92047)

# Data genetation

my_data <- data_generation(n_labeled, N_unlabeled, setting = setting)

Y <- my_data[, 'Y_miss']
S <- my_data[, 'S']

labeled_ind <- which(!is.na(Y))
Y_labeled <- Y[labeled_ind]
S_labeled <- S[labeled_ind]

# Point estimates.

roc_oracle <- supROC(S_labeled = S, Y_labeled = my_data[, "Y"])
roc_sl <- supROC(S_labeled = S_labeled, Y_labeled = Y_labeled)
roc_ss <- ssROC(S, Y)


# Perturbation.

nbt <- 10 # set small for example.
roc_sl_pert <- perturbation(nbt, S, Y, "supROC")
roc_ssl_pert <- perturbation(nbt, S, Y, "ssROC")

# Relative Efficiency

re(roc_sl_pert, roc_ssl_pert, fpr_threshold = 0.1)
```

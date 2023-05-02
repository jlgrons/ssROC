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
setting <- 2

set.seed(92047)

# Data genetation
my_data <- data_generation(n_labeled,
                           N_unlabeled,
                           setting = setting)

Y <- my_data[, 'Y_miss']
S <- my_data[, 'S']

labeled_ind <- which(!is.na(Y))
Y_labeled <- Y[labeled_ind]
S_labeled <- S[labeled_ind]

# Point estimates.

roc_oracle <- supervised(S, my_data[, 'Y'])
roc_sl <- supervised(S_labeled, Y_labeled)
roc_ss <- ssROC(S, Y)


# Perturbation.

nbt <- 2 # set small for example.
roc_sl_pert <- pertubation(nbt, S_labeled, Y_labeled, S, Y, "supervised")
roc_ssl_pert <- pertubation(nbt, S_labeled, Y_labeled, S, Y, "ssROC")
```

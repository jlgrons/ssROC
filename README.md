# ssROC
R package for Semi-Supervised ROC (ssROC) Analysis for Reliable and Streamlined Evaluation of Phenotyping Algorithms.

# Note: This used an old version of the code.

# Installation
```{R, eval = FALSE}
devtools::install_github(repo = "https://github.com/jlgrons/ssROC")
```

# Example
```{R, eval = FALSE}
library(ssROC)

## Set up parameters.
n <- 300
N <- 10000 - 300

set.seed(92047)

## Generate data.
my_data <- data_generation(n, N, setting = 1)

Y <- my_data[, 'Y_miss']
S <- my_data[, 'S']
labeled_ind <- which(!is.na(Y))
Yt <- Y[labeled_ind]
St <- S[labeled_ind]

## Point estimates.
roc.sl0 <- supervised(S, my_data[, 'Y'])
roc.sl <- supervised(St, Yt)
roc.ssROC <- ssROC(S, Y)

## Pertubation resampling for standard error estimation.
roc.sl.pert <- pertubation(nbt = 5, S_labeled = St, Y_labeled = Yt, S = S, Y =Y , method = "supervised")
roc.ssROC.pert <- pertubation(nbt = 5, S_labeled = St, Y_labeled = Yt, S = S, Y = Y, method = "ssROC")
```


devtools::install_github(repo = "https://github.com/ChuanHong/ssROCtesting", force=T)
library(ssROC)
n <- 300
N <- 5000 - 300
p <- 0.3
boot <- FALSE 
nbt <- 2 
m2 <- 0.3
m1 <- 0.7
s2 <- 0.2; s1 <- 0.3
my_data <- data_generation(n, N, m1, s1, m2, s2, p, misspec=F)

Y <- my_data[, 'Y_miss']
S <- my_data[, 'S']
labeled_ind <- which(!is.na(Y))
Yt <- Y[labeled_ind]
St <- S[labeled_ind]

## point estimates
roc.sl0 <- supervised(S, my_data[, 'Y'])
roc.sl <- supervised(St, Yt)
roc.ssROC <- ssROC(S, Y)

## pertubation
roc.sl.pert=pertubation(nbt=nbt, S_labeled=St,Y_labeled=Yt, S=S, Y=Y, method="supervised")
roc.ssROC.pert=pertubation(nbt=nbt, S_labeled=St,Y_labeled=Yt, S=S, Y=Y, method="ssROC")


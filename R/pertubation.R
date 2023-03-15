# -----------------------------------------------------------------------------
# Pertubation bootstrap for standard errors
# -----------------------------------------------------------------------------

#' Pertubation bootstrap for calculating the standard errors of the ROC estimates
#'
#' @param nbt number of replicates for pertubation boostrap
#' @param S_labeled phenotyping score in labeled set
#' @param Y_labeled outcome Y in labeled set
#' @param S phenotyping score in unlabeled set
#' @param Y outcome Y in unlabeled set
#' @param method validation methods, choosing from supervised and ssROC
#' @return List containing AUC and ROC estimates for each perbutation bootstrap
#' @export

pertubation=function(nbt, S_labeled,Y_labeled, S, Y, method){
  res.bt=NULL
  for(ibt in 1:nbt){
    #print(ibt)
    ptb_wgt <- 4*rbeta(length(S_labeled), 1/2, 3/2)
    ptb_wgt_unlabeled <- 4*rbeta(sum(is.na(Y)), 1/2, 3/2)
    
    if(method=="supervised"){res.bt[[ibt]] <- tryCatch(supervised(S_labeled,Y_labeled, W_labeled = ptb_wgt), 
                                                       error=function(e) NA)
    }
    if(method=="ssROC"){res.bt[[ibt]] <- tryCatch(ssROC(S, Y, W_labeled = ptb_wgt,
                                                        W_unlabeled = ptb_wgt_unlabeled),
                                                  error=function(e) NA)
    
    #res.bt[[ibt]]=tryCatch(res.bt[[ibt]][setdiff(ls(res.bt[[ibt]]), "mhat")],error=function(e) NA)
    }
  }
  res.bt
}

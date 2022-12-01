## Confidence Interval and Coverage Probability Functions


logit = function(x){log(x/(1-x))}


normal_based_CI <- function(est, se, level = 0.05){
  upper <- est + qnorm(1-level/2)*se
  lower <- est - qnorm(1-level/2)*se
  return(cbind(lower, upper))
}


logit_trans_normal_based_CI <- function(est, se, level = 0.05){
  upper <- logit(est) + qnorm(1-level/2)*se/(est*(1-est))
  lower <- logit(est) - qnorm(1-level/2)*se/(est*(1-est))
  return(cbind(lower, upper))
}


coverage_probability <- function(lower, upper, truth){
  mean(I(upper >= truth)*I(lower <= truth))
}


get_resamp_sd <- function(num_sims, num_bs_reps, resamps){
  sapply(1:num_sims, function(xx)
    colSds(resamps[(num_bs_reps*(xx-1) + 1):(num_bs_reps*(xx-1) + num_bs_reps), ]))}

## Helper functions for ROC estimation - These need to be cleaned.
         
# -----------------------------------------------------------------------------
# Helper functions for ROC estimation
# -----------------------------------------------------------------------------
#' Helper functions for ROC estimation
#'
#' @param tpr True positive rate
#' @param fpr False positive rate
#' @param p1 Actual postive rate
#' @param S Actual score
#' @param S.ini        
#' @param fpr0 sequennce of fpr        
#'
#' @output
#' @return List containing:
#' \itemize{
#'   \item `auc`, auc
#'   \item `cut`, sequence of cut 
#'   \item `p.pos`, positive rate
#'   \item `fpr`, False positive rate
#'   \item `tpr`, True positive rate
#'   \item `ppv`, Positive predictive value
#'   \item `npv`, Negative predictive value
#'   \item `F.score`, F score test's accuracy        
#' }         
         
         
outpt.FUN = function(tpr, fpr, p1, S, S.ini, fpr0 = seq(0.01, 0.99, by = 0.01)){
  cuts=unique(sort(S))
  p0=1-p1
  ppv = tpr*p1/(tpr*p1+fpr*p0)
  npv = (1-fpr)*p0/((1-fpr)*p0+(1-tpr)*p1)
  fscore = 2*ppv*tpr/(ppv+tpr)
  p.pos = unlist(lapply(1:length(cuts), function(ll) mean(S>=cuts[ll])))
  cuts = quantile(S.ini, 1-p.pos)
  junk = cbind("cut"= cuts, "p.pos"=p.pos, "fpr"=round(fpr, 2),
               "tpr"=tpr, "ppv"=ppv, "npv"=npv, "F.score"=fscore)
  id.print = unlist(lapply(fpr0, function(x) which.min(abs(x-fpr))[1]))
  auc = auc.FUN(tpr, fpr)
  out = list(auc = auc, roc = junk[id.print,])
  out
}

#' inverse logit 
g.logit = function(xx){exp(xx)/(exp(xx)+1)}

# -----------------------------------------------------------------------------
# Create VTM
# -----------------------------------------------------------------------------
#' Creat VTM
                           
#' @ param vc vector
#' @ param dm dimension
#'
#' @ return vector transmission matrix                          
                           
VTM<-function(vc, dm){
  matrix(vc, ncol=length(vc), nrow=dm, byrow=T)
}

# -----------------------------------------------------------------------------
# 
# -----------------------------------------------------------------------------
#' 
#'
#' @param yy          
#' @param Yi Observed outcome
#' @param Di Predicted outcome
#' @param yes,smooth                          
                                                     
S.FUN <- function(yy,Yi,Di,yes.smooth=F)
{
  if(yes.smooth){
    Y1i = Yi[Di==1]; n1 = sum(Di); bw = bw.nrd(Y1i)/n1^0.6
    c(t(rep(1/n1,n1))%*%pnorm((Y1i-VTM(yy,n1))/bw))
  }else{
    return((sum.I(yy,"<",Yi,Vi=Di)+sum.I(yy,"<=",Yi,Vi=Di))/sum(Di)/2)
  }
}


# -----------------------------------------------------------------------------
# 
# -----------------------------------------------------------------------------
#' 
#' @param yy 
#' @param FUN relation between yy and Yi
#' @param Yi Observed outcome
#' @param Vi target vector
                           
sum.I <- function(yy,FUN,Yi,Vi=NULL)
  # Vi weight
{
  if (FUN=="<"|FUN==">=") { yy <- -yy; Yi <- -Yi}
  # for each distinct ordered failure time t[j], number of Xi < t[j]
  pos <- rank(c(yy,Yi),ties.method='f')[1:length(yy)]-rank(yy,ties.method='f')
  if (substring(FUN,2,2)=="=") pos <- length(Yi)-pos # number of Xi>= t[j]
  if (!is.null(Vi)) {
    ## if FUN contains '=', tmpind is the order of descending
    if(substring(FUN,2,2)=="=") tmpind <- order(-Yi) else  tmpind <- order(Yi)
    Vi <- apply(as.matrix(Vi)[tmpind,,drop=F],2,cumsum)
    return(rbind(0,Vi)[pos+1,])
  } else return(pos)
}

# -----------------------------------------------------------------------------
# 
# -----------------------------------------------------------------------------
#' 
#' @param uu set of numerica value speckfying wgere interpolation take place
#' @param Yi Observed outcome
#' @param Di Predicted oputcome
#' @param yes.smooth 
#'                           
#' @return list y

Sinv.FUN <- function(uu,Yi,Di,yes.smooth=F)
{
  yy0<-unique(sort(Yi,decreasing=T)); ss0 <- S.FUN(yy0,Yi,Di,yes.smooth=yes.smooth)
  return(approx(ss0[!duplicated(ss0)],yy0[!duplicated(ss0)],uu,method="linear",f=0,rule=2)$y)
}

# -----------------------------------------------------------------------------
# AUC function
# -----------------------------------------------------------------------------
#' AUC function
#'
#' @param tpr True positive rate
#' @param fpr False positive rate  
#'
#' @return AUC                           

auc.FUN=function(TPR, FPR){
  sens=c(sort(TPR,decreasing=T),0); omspec=c(sort(FPR, decreasing=T),0)
  height = (sens[-1]+sens[-length(sens)])/2
  width = -diff(omspec)
  AUC = sum(height*width)
  AUC
}

# -----------------------------------------------------------------------------
# 
# -----------------------------------------------------------------------------
#'
#'
#' @param St 
#' @param Yt
#' @param Sv 
#' @param bw                        
#' @param Wt 
#' @param kern.mat
                           
NP.REG <- function(St, Yt, Sv, bw, Wt = NULL, kern.mat=NULL){
  nv <- length(Sv)
  nt <- length(St)
  if (is.null(Wt)){Wt <- rep(1,nt)}
  if(is.null(kern.mat)){kern.mat <- sapply(1:nt, function(kk) dnorm(Sv - rep(St[kk], nv), sd = bw))}
  nw.est =  kern.mat %*% (Yt*Wt) * (1/(kern.mat%*%Wt));
  return(nw.est);
}


## Helper functions for cutoff estimation -- need to be cleaned


# -----------------------------------------------------------------------------
# choosing cut
# -----------------------------------------------------------------------------
#'choosing cut   
#'                                           
#' @param roc roc table
#' @param col.nm quantity of interest char
#' @param value desired value                                        
                                           
newfunction=function(roc,col.nm,value){
  cutoff.line = roc[which.min(abs(roc[,col.nm]-value)), ]
  cutoff = cutoff[cut]
}

# -----------------------------------------------------------------------------
# choosing cut
# -----------------------------------------------------------------------------
#'choosing cut                                            

#' @param roc : num [1:99,1:6] roc table
#' @param col.nm : quantity of interest char
#' @param value : real number
                                           
cutoff.choose=function(roc,col.nm,value){
  roc = data.frame(roc)
  cutoff.line = roc[which.min(abs(roc[,col.nm]-value))[1], ]
  cutoff = cutoff.line$cut
  res = list(cutoff,cutoff.line)
  res
}

#Interpolated ROC function used by supervised and ss imp estimators.
# -----------------------------------------------------------------------------
# Interpolated ROC function used by supervised and ss imp estimators.
# -----------------------------------------------------------------------------
#'Interpolated ROC function used by supervised and ss imp estimators.
#'                                           
#' @param S Score
#' @param Y Outcome
#' @param w weight vector
#' @param fpr_vals equennce of fpr
#'                                           
#' @return List containing:
#' \itemize{
#'   \item `cut`, sequence of cut 
#'   \item `p.pos`, positive rate
#'   \item `fpr`, False positive rate
#'   \item `tpr`, True positive rate
#'   \item `ppv`, Positive predictive value
#'   \item `npv`, Negative predictive value
#'   \item `F.score`, F score test's accuracy
#'   \item `auc`, auc                                           
#' }                                       
                                           
interpolated_ROC <- function(S, Y, W = NULL,
                             fpr_vals = seq(0.01, 0.99, by = 0.01)){
  
  cuts <- unique(sort(S))
  nc <- length(cuts)
  
  p_pos <- unlist(lapply(1:length(cuts), function(ll) 
    sum((S>=cuts[ll]) * W) / sum(W))
  )
  
  TPR_c <- sum.I(cuts, "<=", S, Y * W) / sum(Y * W)
  FPR_c <- sum.I(cuts, "<=", S, (1-Y) * W) / sum((1-Y) * W)
  
  auc <- sum(TPR_c[-1] * (FPR_c[-nc]-FPR_c[-1]))
  
  mu1 <- sum(Y * W) / sum(W)
  mu0 <- 1 - mu1
  
  PPV_c <- (TPR_c * mu1) / (TPR_c * mu1 + FPR_c * mu0)
  NPV_c <- ((1-FPR_c) * mu0) / ((1-FPR_c) * mu0 + (1-TPR_c) * mu1)
  
  fscore_c <- (2 * PPV_c * TPR_c) / (PPV_c + TPR_c)
  
  roc_c <- cbind("cut"= cuts,"p.pos" = p_pos, "FPR"=FPR_c,"TPR"=TPR_c,
                 "PPV"=PPV_c,"NPV"=NPV_c, "F.score" = fscore_c)
  
  roc <- sapply(1:ncol(roc_c), function(kk){approx(roc_c[,"FPR"], roc_c[,kk],
                                                   fpr_vals, rule = 2)$y})
  colnames(roc) <-colnames(roc_c)
  
  return(list(roc = roc, auc = auc))
}






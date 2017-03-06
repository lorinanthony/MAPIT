#This function will run a version of MAPIT under the following variations:
#(1) Standard Model ---> y = m+g+e
#(2) Standard + Covariate Model where Z is a matrix of covariates ---> y = Zu+m+g+e
#(3) Standard + Common Environment Model where c controls for extra environmental effects and population structure ---> y = m+g+c+e
#(3) Standard + Covariate + Common Environment Model ---> y = Zu+m+g+c+e
#
#And will consider the following two hypothesis testing strategies:
#(1) The Normal or Z-test
#(2) Davies Method
#
#Default: Standard Model under the Normal or Z-test
#
#In order to compile, this wrapper function requires the following R libraries:
library(doParallel)
library(Rcpp)
library(RcppArmadillo)
library(CompQuadForm)

######################################################################################

MAPIT = function(X, y, Z = NULL,C = NULL,threshold = NULL){
  ### Run the analysis using the Normal/Z-test ###
  
  ### Run Model 1 by Default ###
  if(is.null(Z)==TRUE&&is.null(C)==TRUE){
    vc.mod = MAPIT1_Normal(X,y)
  }
  ### Test to Run Model 2 ###
  if(is.null(Z)==FALSE&&is.null(C)==TRUE){
    vc.mod = MAPIT2_Normal(X,y,Z)
  }
  ### Test to Run Model 3 ###
  if(is.null(Z)==TRUE&&is.null(C)==FALSE){
    vc.mod = MAPIT3_Normal(X,y,C)
  }
  ### Test to Run Model 4 ###
  if(is.null(Z)==FALSE&&is.null(C)==FALSE){
    vc.mod = MAPIT4_Normal(X,y,Z,C)
  }
  
  ### Record the Results of the Normal Z-Test ###
  pvals = vc.mod$pvalues; names(pvals) = rownames(X)
  pves = vc.mod$PVE; names(pves) = rownames(X)
  
  ### Re-run the analysis under the Davies Method ###
  if(is.null(threshold)==FALSE){
    
    ### Find the indices of the p-values that are below 0.05 ###
    ind = which(pvals<=threshold)
    
    ### Run Model 1 by Default ###
    if(is.null(Z)==TRUE&&is.null(C)==TRUE){
      vc.mod = MAPIT1_Davies(X,y,ind)
    }
    ### Test to Run Model 2 ###
    if(is.null(Z)==FALSE&&is.null(C)==TRUE){
      vc.mod = MAPIT2_Davies(X,y,Z,ind)
    }
    ### Test to Run Model 2 ###
    if(is.null(Z)==TRUE&&is.null(C)==FALSE){
      vc.mod = MAPIT3_Davies(X,y,C,ind)
    }
    ### Test to Run Model 2 ###
    if(is.null(Z)==FALSE&&is.null(C)==FALSE){
      vc.mod = MAPIT4_Davies(X,y,Z,C,ind)
    }
    ### Apply Davies Exact Method ###
    vc.ts = vc.mod$Est
    names(vc.ts) = rownames(X)
    
    davies.pvals = c()
    for(i in 1:length(vc.ts)){
      lambda = sort(vc.mod$Eigenvalues[,i],decreasing = T)
      
      Davies_Method = davies(vc.mod$Est[i], lambda = lambda, acc=1e-8)
      davies.pvals[i] = 2*min(Davies_Method$Qq, 1-Davies_Method$Qq)
      names(davies.pvals)[i] = names(vc.ts[i])
    }
    pvals[ind] = davies.pvals[ind]
  }
  return(list("pvalues"=pvals,"pves"=pves))
}
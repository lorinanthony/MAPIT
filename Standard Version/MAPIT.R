###***Parallelized Version of the MArginal ePIstasis Test (MAPIT)***###

#This function will run a version of the MArginal ePIstasis Test (MAPIT) under the following model variations:

#(1) Standard Model ---> y = m+g+e where m ~ MVN(0,omega^2K), g ~ MVN(0,sigma^2G), e ~ MVN(0,tau^2M). Recall from Crawford et al. (2017) that m is the combined additive effects from all other variants, and effectively rep- resents the additive effect of the kth variant under the polygenic background of all other variants; K = X_{−k}^tX_{−k}/(p − 1) is the genetic relatedness matrix computed using genotypes from all variants other than the kth; g is the summation of all pairwise interaction effects between the kth variant and all other variants; G = DKD represents a relatedness matrix computed based on pairwise interaction terms between the kth variant and all other variants. Here, we also denote D = diag(x_k) to be an n × n diagonal matrix with the genotype vector x_k as its diagonal elements. It is important to note that both K and G change with every new marker k that is considered. Lastly; M is a variant specific projection matrix onto both the null space of the intercept and the corresponding genotypic vector x_k.

#(2) Standard + Covariate Model ---> y = Wa+m+g+e where W is a matrix of covariates with effect sizes a.

#(3) Standard + Common Environment Model ---> y = m+g+c+e where c ~ MVN(0,eta^2C) controls for extra environmental effects and population structure with covariance matrix C.

#(4) Standard + Covariate + Common Environment Model ---> y = Wa+m+g+c+e

#This function will consider the following three hypothesis testing strategies which are featured in Crawford et al. (2017):
#(1) The Normal or Z-test
#(2) Davies Method
#(3) Hybrid Method (Z-test + Davies Method)

### Review of the function parameters ###
#'X' is the pxn genotype matrix where p is the number of variants and n is the number of samples. Must be a matrix and not a data.frame.
#'y' is the nx1 vector of quantitative or continuous traits.
#'W' is the matrix qxn matrix of covariates. Must be a matrix and not a data.frame.
#'C' is an nxn covariance matrix detailing environmental effects and population structure effects.
#'hybrid' is a parameter detailing if the function should run the hybrid hypothesis testing procedure between the normal Z-test and the Davies method. Default is TRUE.
#'threshold' is a parameter detailing the value at which to recalibrate the Z-test p-values. If nothing is defined by the user, the default value will be 0.05 as recommended by the Crawford et al. (2017). 
#'test' is a parameter defining what hypothesis test should be implemented. Takes on values 'normal' or 'davies'. This parameter only matters when hybrid = FALSE. If test is not defined when hybrid = FALSE, the function will automatically use test = 'normal'.
#'k' is a parameter specifying the known prevelance of a disease or trait in the population as in Crawford and Zhou (2018). This parameter is only used if y is a binary class of labels (0,1). If this parameter is not set then the default k = 0.5 will used.

######################################################################################
######################################################################################
######################################################################################

MAPIT = function(X, y, W = NULL,C = NULL,hybrid = TRUE,threshold = 0.05,test = NULL,k = NULL){
  ### Install the necessary libraries ###
  usePackage("doParallel")
  usePackage("Rcpp")
  usePackage("RcppArmadillo")
  usePackage("CompQuadForm")
  usePackage("truncnorm")
  
  ### Check to See if we should run LT-MAPIT ###
  if(sum(y%in%c(0,1))==length(y)){
    if(k = NULL){warning("The liability threshold model is going to be used but no disease prevelance is defined! Using the default 0.5.");k = 0.5}
    ### Save the labels ###
    cc = y
    
    ### Find the Number of Cases and Controls ###
    n.cases = sum(cc==1)
    n.controls = sum(cc==0)
    
    ### Set the Threshold ###
    thresh=qnorm(1-k,mean=0,sd=1)
    
    ### Bernoulli Distributed Case and Control Data ###
    Pheno=rep(NA,length(cc));
    Pheno[cc==0] = rtruncnorm(n.controls,b=thresh)
    Pheno[cc==1] = rtruncnorm(n.cases,a=thresh)
    y = Pheno
  }
  
  ### Run the analysis using the hybrid Test ###
  if(hybrid==TRUE){
    ### Run Model 1 by Default ###
    if(is.null(W)==TRUE&&is.null(C)==TRUE){
      vc.mod = MAPIT1_Normal(X,y)
    }
    ### Test to Run Model 2 ###
    if(is.null(W)==FALSE&&is.null(C)==TRUE){
      vc.mod = MAPIT2_Normal(X,y,W)
    }
    ### Test to Run Model 3 ###
    if(is.null(W)==TRUE&&is.null(C)==FALSE){
      vc.mod = MAPIT3_Normal(X,y,C)
    }
    ### Test to Run Model 4 ###
    if(is.null(W)==FALSE&&is.null(C)==FALSE){
      vc.mod = MAPIT4_Normal(X,y,W,C)
    }
    ### Record the Results of the Normal Z-Test ###
    pvals = vc.mod$pvalues; names(pvals) = rownames(X)
    pves = vc.mod$PVE; names(pves) = rownames(X)
    
    ### Find the indices of the p-values that are below 0.05 ###
    ind = which(pvals<=threshold)
    
    ### Run Model 1 by Default ###
    if(is.null(W)==TRUE&&is.null(C)==TRUE){
      vc.mod = MAPIT1_Davies(X,y,ind)
    }
    ### Test to Run Model 2 ###
    if(is.null(W)==FALSE&&is.null(C)==TRUE){
      vc.mod = MAPIT2_Davies(X,y,W,ind)
    }
    ### Test to Run Model 2 ###
    if(is.null(W)==TRUE&&is.null(C)==FALSE){
      vc.mod = MAPIT3_Davies(X,y,C,ind)
    }
    ### Test to Run Model 2 ###
    if(is.null(W)==FALSE&&is.null(C)==FALSE){
      vc.mod = MAPIT4_Davies(X,y,W,C,ind)
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
  
  ### Run the analysis using the either normal test or Davies method ###
  if(hybrid==FALSE&&is.null(test)==TRUE){
    warning("No test specified for hypothesis testing. Using normal Z-test (default).")
    ### Run Model 1 by Default ###
    if(is.null(W)==TRUE&&is.null(C)==TRUE){
      vc.mod = MAPIT1_Normal(X,y)
    }
    ### Test to Run Model 2 ###
    if(is.null(W)==FALSE&&is.null(C)==TRUE){
      vc.mod = MAPIT2_Normal(X,y,W)
    }
    ### Test to Run Model 3 ###
    if(is.null(W)==TRUE&&is.null(C)==FALSE){
      vc.mod = MAPIT3_Normal(X,y,C)
    }
    ### Test to Run Model 4 ###
    if(is.null(W)==FALSE&&is.null(C)==FALSE){
      vc.mod = MAPIT4_Normal(X,y,W,C)
    }
    
    ### Record the Results of the Normal Z-Test ###
    pvals = vc.mod$pvalues; names(pvals) = rownames(X)
    pves = vc.mod$PVE; names(pves) = rownames(X)
  }
  
  ### Run the analysis using the either normal test or Davies method ###
  if(hybrid==FALSE&&is.null(test)==FALSE){
    if(test=="normal"||test=="Normal"){
      ### Run Model 1 by Default ###
      if(is.null(W)==TRUE&&is.null(C)==TRUE){
        vc.mod = MAPIT1_Normal(X,y)
      }
      ### Test to Run Model 2 ###
      if(is.null(W)==FALSE&&is.null(C)==TRUE){
        vc.mod = MAPIT2_Normal(X,y,W)
      }
      ### Test to Run Model 3 ###
      if(is.null(W)==TRUE&&is.null(C)==FALSE){
        vc.mod = MAPIT3_Normal(X,y,C)
      }
      ### Test to Run Model 4 ###
      if(is.null(W)==FALSE&&is.null(C)==FALSE){
        vc.mod = MAPIT4_Normal(X,y,W,C)
      }
      
      ### Record the Results of the Normal Z-Test ###
      pvals = vc.mod$pvalues; names(pvals) = rownames(X)
      pves = vc.mod$PVE; names(pves) = rownames(X)
    }
    
    if(test=="davies"||test=="Davies"){
      ### Find the indices of the p-values that are below 0.05 ###
      ind = 1:nrow(X)
      
      ### Run Model 1 by Default ###
      if(is.null(W)==TRUE&&is.null(C)==TRUE){
        vc.mod = MAPIT1_Davies(X,y,ind)
      }
      ### Test to Run Model 2 ###
      if(is.null(W)==FALSE&&is.null(C)==TRUE){
        vc.mod = MAPIT2_Davies(X,y,W,ind)
      }
      ### Test to Run Model 2 ###
      if(is.null(W)==TRUE&&is.null(C)==FALSE){
        vc.mod = MAPIT3_Davies(X,y,C,ind)
      }
      ### Test to Run Model 2 ###
      if(is.null(W)==FALSE&&is.null(C)==FALSE){
        vc.mod = MAPIT4_Davies(X,y,W,C,ind)
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
      pvals = davies.pvals; pves = vc.mod$PVE
    }
  }
  
  
  return(list("pvalues"=pvals,"pves"=pves))
}

usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE)
  require(p, character.only = TRUE)
}
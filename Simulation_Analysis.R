### Clear Console ###
cat("\014")

### Clear Environment ### 
rm(list = ls(all = TRUE))

### Load in the R libraries ###
library(coda)
library(MASS)
library(doParallel)
library(Rcpp)
library(RcppArmadillo)
library(RcppParallel)
library(BGLR)
library(monomvn)
library(kernlab)
library(CompQuadForm)

### Load in functions to make QQ-plot plots ###
source("QQplot.R")

######################################################################################
######################################################################################
######################################################################################

### Set the Seed for the analysis ###
#set.seed(11151990)

# Data are p single nucleotide polymorphisms (SNPs) with simulated genotypes.
# Simulation Parameters: 
# (1) ind = # of samples 
# (2) nsnp = number of SNPs or variants
# (3) PVE = phenotypic variance explained/broad-sense heritability (H^2)
# (4) rho = measures the portion of H^2 that is contributed by the marignal (additive) effects

ind = 3e3; nsnp = 1e4; pve=0.6; rho=0.5;

### Simulate the genotypes such that all variants have minor allele frequency (MAF) > 0.05 ###
# NOTE: As in the paper, we center and scale each genotypic vector such that every SNP has mean 0 and standard deviation 1.
maf <- 0.05 + 0.45*runif(nsnp)
X   <- (runif(ind*nsnp) < maf) + (runif(ind*nsnp) < maf)
X   <- matrix(as.double(X),ind,nsnp,byrow = TRUE)
Xmean=apply(X, 2, mean); Xsd=apply(X, 2, sd); X=t((t(X)-Xmean)/Xsd)

### Set the size of causal SNP groups ###
ncausal1= 10 # Set 1 of causal SNPs 
ncausal2 = 10 # Set 2 of Causal SNPs
ncausal3 = 1e3-ncausal1-ncausal2 # Set 3 of Causal SNPs with only marginal effects

# NOTE: As in the paper, we randomly choose causal variants that classify into three groups: (i) a small set of interaction SNPs, (ii) a larger set of interaction SNPs, and (iii) a large set of additive SNPs. In the simulations carried out in this study, SNPs interact between sets, so that SNPs in the first group interact with SNPs in the second group, but do not interact with variants in their own group (the same applies to the second group). One may view the SNPs in the first set as the “hubs” in an interaction map. We are reminded that interaction (epistatic) effects are different from additive effects. All causal SNPs in both the first and second groups have additive effects and are involved in pairwise interactions, while causal SNPs in the third set only have additive effects.

### Select causal SNPs in groups 1 and 2 ###
snp.ids = 1:nsnp
s1=sample(snp.ids, ncausal1, replace=F)
s2=sample(snp.ids[-s1], ncausal2, replace=F)
s3=sample(snp.ids[-c(s1,s2)], ncausal3, replace=F)

Xcausal1=X[,s1]; Xcausal2=X[,s2]; Xcausal3=X[,s3]

### Create the epistatic design matrix ###
W=c()
for(i in 1:ncausal1){
  W=cbind(W,Xcausal1[,i]*Xcausal2)
}
dim(W)

### Simulate marginal (additive) effects ###
Xmarginal=cbind(Xcausal1,Xcausal2,Xcausal3)
beta=rnorm(dim(Xmarginal)[2])
y_marginal=Xmarginal%*%beta
beta=beta*sqrt(pve*rho/var(y_marginal))
y_marginal=Xmarginal%*%beta

### Simulate epistatic effects ###
alpha=rnorm(dim(Xepi)[2])
y_epi=W%*%alpha
alpha=alpha*sqrt(pve*(1-rho)/var(y_epi))
y_epi=W%*%alpha

# NOTE: Each effect size for the marginal, epistatic, and random error effects are drawn from a standard normal distribution. Meaning beta ~ MVN(0,I), alpha ~ MVN(0,I), and epsilon ~ MVN(0,I). We then scale both the additive and pairwise genetic effects so that collectively they explain a fixed proportion of genetic variance. Namely, the additive effects make up rho%, while the pairwise interactions make up the remaining (1 − rho)%. Once we obtain the final effect sizes for all causal SNPs, we draw errors to achieve the target H^2.

### Simulate residual error ###
y_err=rnorm(ind)
y_err=y_err*sqrt((1-pve)/var(y_err))

# Simulate the continuous phenotypes
y = y_marginal+y_epi+y_err

### Check dimensions and add SNP names ###
dim(X); dim(y)
colnames(X) = paste("SNP",1:ncol(X),sep="")

### Save the names of the causal SNPs ###
SNPs = colnames(X)[c(s1,s2)]

######################################################################################
######################################################################################
######################################################################################

### Running MEPIT using the Normal Method ###

### Load in the C++ functions ###
sourceCpp("MEPIT_Functions.cpp")

Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
cores = detectCores() 
ptm <- proc.time() #Start clock
vc.mod = MEPIT_Normal(t(X),y,cores=cores)
proc.time() - ptm #Stop clock

vc.pvals1 = vc.mod$pvalues
names(vc.pvals1) = colnames(X)

######################################################################################
######################################################################################
######################################################################################

### Find the indices of the p-values that are below 0.05 ###
ind = which(vc.pvals1<= 0.05)

### Rerun MEPIT using the Davies Method to Recalibrate the p-values for SNPs <= 0.05 ###
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
cores = detectCores() 
ptm <- proc.time() #Start clock
vc.mod = MEPIT_Davies(t(X),y,ind,cores=cores)
proc.time() - ptm #Stop clock

### Apply Davies Exact Method ###
vc.ts = vc.mod$Est
names(vc.ts) = colnames(X)

vc.pvals2 = c()
for(i in 1:length(vc.ts)){
  lambda = sort(vc.mod$Eigenvalues[,i],decreasing = T)
  
  Davies_Method = davies(vc.mod$Est[i], lambda = lambda, acc=1e-8)
  vc.pvals2[i] = 2*min(Davies_Method$Qq, 1-Davies_Method$Qq)
  names(vc.pvals2)[i] = names(vc.ts[i])
}

### Combine the two sets of p-values
vc.pvals = vc.pvals1; vc.pvals[ind] = vc.pvals2[ind]

### Plot observed the p-values on a QQ-plot ###
ggd.qqplot(vc.pvals)
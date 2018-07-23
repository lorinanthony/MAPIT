### Illustrating the Liability Threshold MArginal ePIstasis Test (LT-MAPIT) with Simulations ###

### Clear Console ###
cat("\014")

### Clear Environment ### 
rm(list = ls(all = TRUE))

### Load in the R libraries ###
library(doParallel)
library(Rcpp)
library(RcppArmadillo)
library(RcppParallel)
library(CompQuadForm)
library(truncnorm)

### Load in functions to make QQ-plot plots ###
source("QQplot.R")

#NOTE: This code assumes that the basic C++ functions are set up on the computer in use. If not, the MAPIT functions and Rcpp packages will not work properly. Mac users please refer to the homebrew applications and install the gcc commands listed in the README.md file before running the rest of the code [Warning: This step may take about an hour...].

### Load in the C++ MAPIT wrapper functions ###
source("MAPIT.R"); sourceCpp("MAPIT.cpp")

######################################################################################
######################################################################################
######################################################################################

### Set the Seed for the analysis ###
#set.seed(11151990)

# Data are p single nucleotide polymorphisms (SNPs) with simulated genotypes.
# Simulation Parameters: 
# (1) ind = # of samples in the population
# (2) nsnp = number of SNPs or variants
# (3) PVE = phenotypic variance explained/broad-sense heritability (H^2)
# (4) rho = measures the portion of H^2 that is contributed by the marignal (additive) effects
# (5) k = assumed disease prevelance in the population
# (6) samp = # of samples to analyze

ind = 1e6; nsnp = 1e4; pve=0.6; rho=0.5; k = 0.001; samp = 500

### Simulate the genotypes such that all variants have minor allele frequency (MAF) > 0.05 ###
# NOTE: As in the paper, we center and scale each genotypic vector such that every SNP has mean 0 and standard deviation 1.
maf <- 0.05 + 0.45*runif(nsnp)
Geno   <- (runif(ind*nsnp) < maf) + (runif(ind*nsnp) < maf)
Geno   <- matrix(as.double(Geno),ind,nsnp,byrow = TRUE)
Xmean=apply(Geno, 2, mean); Xsd=apply(Geno, 2, sd); Geno=t((t(Geno)-Xmean)/Xsd)

### Set the size of causal SNP groups ###
ncausal1= 10 # Set 1 of causal SNPs 
ncausal2 = 10 # Set 2 of Causal SNPs
ncausal3 = 1e3-ncausal1-ncausal2 # Set 3 of Causal SNPs with only marginal effects

# NOTE: As in the paper, we randomly choose causal variants that classify into three groups: (i) a small set of interaction SNPs, (ii) a larger set of interaction SNPs, and (iii) a large set of additive SNPs. In the simulations carried out in this study, SNPs interact between sets, so that SNPs in the first group interact with SNPs in the second group, but do not interact with variants in their own group (the same applies to the second group). One may view the SNPs in the first set as the “hubs” in an interaction map. We are reminded that interaction (epistatic) effects are different from additive effects. All causal SNPs in both the first and second groups have additive effects and are involved in pairwise interactions, while causal SNPs in the third set only have additive effects.

#Select Causal SNPs
s = 1:nsnp
s1=sample(s, ncausal1, replace=F)
s2=sample(s[s%in%s1==FALSE], ncausal2, replace=F)
s3=sample(s[s%in%c(s1,s2)==FALSE], ncausal3, replace=F)

#Create Causal Epistatic Matrix
Xcausal1=Geno[,s1]; Xcausal2=Geno[,s2]; Xcausal3=Geno[,s3]
Xepi=c()
for(i in 1:ncausal1){
  Xepi=cbind(Xepi,Xcausal1[,i]*Xcausal2)
}
dim(Xepi)

#Marginal effects only
Xmarginal=cbind(Xcausal1,Xcausal2,Xcausal3)
beta=rnorm(dim(Xmarginal)[2])
z_marginal=c(Xmarginal%*%beta)
beta=beta*sqrt(pve*rho/var(z_marginal))
z_marginal=Xmarginal%*%beta

#Epistatic effects
beta=rnorm(dim(Xepi)[2])
z_epi=c(Xepi%*%beta)
beta=beta*sqrt(pve*(1-rho)/var(z_epi))
z_epi=Xepi%*%beta

# NOTE: Each effect size for the marginal, epistatic, and random error effects are drawn from a standard normal distribution. Meaning beta ~ MVN(0,I), alpha ~ MVN(0,I), and epsilon ~ MVN(0,I). We then scale both the additive and pairwise genetic effects so that collectively they explain a fixed proportion of genetic variance. Namely, the additive effects make up rho%, while the pairwise interactions make up the remaining (1 − rho)%. Once we obtain the final effect sizes for all causal SNPs, we draw errors to achieve the target H^2.

#Error
z_err=rnorm(ind)
z_err=z_err*sqrt((1-pve)/var(z_err))

#Latent Variables
z=z_marginal+z_epi+z_err

# NOTE: Now that we have the liabilities, we can assign case-control labels according to the disease prevelance parameter. We will treat this like the LT-MAPIT paper and take an equal number of cases and controls.

### Set the Threshold ###
thresh=qnorm(1-k,mean=0,sd=1)

### Find the Number of Cases and Controls ###
n.cases = sum(z>thresh); n.cases/length(z)
n.controls = sum(z<=thresh); n.controls/length(z)

#Bernoulli Distributed Case and Control Data
Pheno=rep(NA,ind); 
Pheno[z<=thresh] = rtruncnorm(n.controls,b=thresh) 
Pheno[z>thresh] = rtruncnorm(n.cases,a=thresh)

### Subsample a particular number of cases and controls ###
cases = sample(which(z>thresh),samp/2,replace = FALSE)
controls = sample(which(z<=thresh),samp/2,replace = FALSE)
y = Pheno[c(cases,controls)]
X = Geno[c(cases,controls),]

### Remove the variables that are not being used ###
rm(z); rm(z_marginal); rm(z_epi); rm(Xmarginal); rm(Xepi)
rm(Xcausal1); rm(Xcausal2); rm(Xcausal3)

### Check dimensions and add SNP names ###
dim(X); length(y)
colnames(X) = paste("SNP",1:ncol(X),sep="")

### Save the names of the causal SNPs ###
SNPs = colnames(X)[c(s1,s2)]

######################################################################################
######################################################################################
######################################################################################

### Running LT-MAPIT using Davies Method ###

### Set the number of cores ###
cores = detectCores()

### Run LT-MAPIT ###
ptm <- proc.time() #Start clock
mapit = MAPIT(t(X),y,hybrid=FALSE,test="davies",cores=cores)
proc.time() - ptm #Stop clock

davies.pvals = mapit$pvalues
names(davies.pvals) = colnames(X)

### NOTE: The first set of p-values will be deflated (as mentioned in Crawford and Zhou (2018)). To fix this issues we next run the test statistic readjustment procedure ###
pvals = davies.pvals[davies.pvals>=0]; summary(pvals)
q = qchisq(pvals,df=1,lower.tail = FALSE); summary(q)
gs = median(q); m = 0.4549; summary((q/gs)*m)
pvals = pchisq((q/gs)*m,df=1,lower.tail = FALSE)
summary(pvals)

### Plot observed the p-values on a QQ-plot ###
ggd.qqplot(pvals)

### Look at the causal SNPs in group 1 and 2 ###
pvals[s1]; pvals[s2]

######################################################################################
######################################################################################
######################################################################################

### Running an Informed Exhaustive Search ###

#NOTE: Now we may take only the significant SNPs according to their marginal epistatic effects and run a simple exhaustive search between them

thresh = 0.05/length(pvals) #Set a significance threshold
v = pvals[pvals<=thresh] #Call only marginally significant SNPs
pwise = c()
for(k in 1:length(v)){
  fit = c(); m = 1:length(v)
  for(w in m[-k]){
    pt = lm(y~X[,names(v)[k]]:X[,names(v)[w]])
    fit[w] = coefficients(summary(pt))[8]
    names(fit)[w] = paste(names(v)[k],names(v)[w],sep = "-")
  }
  pwise = c(pwise,fit)
}

### Look at the exhaustive search results ###
vc.pwise = sort(pwise[!is.na(pwise)]) #Get rid of the NAs
vc.pwise = vc.pwise[seq(1,length(vc.pwise),2)] #Only keep the unique pairs
sort(vc.pwise)[1:150] #Sort the pairs in order of significance
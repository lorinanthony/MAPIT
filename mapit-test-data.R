# This file needs to be run in the directory of the
# original MAPIT implementation.

library(doParallel)
library(Rcpp)
library(RcppArmadillo)
library(RcppParallel)
library(CompQuadForm)
library(Matrix)
library(MASS)
library(truncnorm)

# Load in the MAPIT Functions
# sourceCpp("~/git/MAPIT/Standard Version/MAPIT.cpp")
# source("~/git/MAPIT/Standard Version/MAPIT.R")
sourceCpp("~/git/MAPIT/OpenMP Version/MAPIT_OpenMP.cpp")
source("~/git/MAPIT/OpenMP Version/MAPIT_OpenMP.R")

# Set the OpenMP Envrionment
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")

# test <- c("normal", "davies")

# Read in the Data
Data = readRDS("/users/jstamp1/simulations/data/control/null_H6_R10_P10_N1000.rds")
sample <- Data[[1]]
X = sample$genotype
X <- X[,which(apply(X, 2, var) != 0)]
X <- X[, which(apply(X, 2, mean) >= 0.05)]
Xmean=apply(X, 2, mean); Xsd=apply(X, 2, sd); X=t((t(X)-Xmean)/Xsd)

Y = sample$phenotype
y <- Y[,1]

# Set the number of cores
cores = detectCores()

# Run MAPIT: Full Version with Satterthwaite Method
ptm <- proc.time() #Start clock
vc.mod.normal = MAPIT(t(X),y,hybrid=FALSE,test="normal",cores = cores)
vc.mod.davies = MAPIT(t(X),y,hybrid=FALSE,test="davies",cores = cores)
# vc.mod = MAPIT(t(X),y,hybrid=FALSE,test=test)
proc.time() - ptm #Stop clock

results <- list("genotype" = X,
                "phenotype" = y,
                "normal" = vc.mod.normal,
                "davies" = vc.mod.davies)
saveRDS(results, "original_MAPIT.rds")


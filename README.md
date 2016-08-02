# The Marginal EPIstasis Test (MEPIT)

### Introduction
Epistasis, commonly defined as the interaction between multiple genes, is an important genetic component underlying phenotypic variation. Many statistical methods have been developed to model and identify epistatic interactions between genetic variants. However, because of the large combinatorial search space of interactions, most epistasis mapping methods face enormous computational challenges and often suffer from low statistical power. In [Crawford, Mukherjee, and Zhou (2016)](http://biorxiv.org/content/early/2016/07/31/066985), we present a novel, alternative strategy for mapping epistasis: the marginal EPIstasis test (MEPIT). Our method examines one variant at a time, and estimates and tests its "marginal epistatic effects" --- the combined pairwise interaction effects between a given variant and all other variants. By avoiding explicitly searching for interactions, our method avoids the large combinatorial search space and improves power. Our method is novel and relies on a recently developed variance component estimation method for efficient and robust parameter inference and p-value computation.

MEPIT is implemented as a set of R and C++ routines, which can be carried out within an R environment.


### The R Environment
R is a widely-used, free and open source software environment for statistical computing and graphics. The most recent version of R can be downloaded from the 
[Comprehensive R Archive Network (CRAN)](http://cran.r-project.org/)
CRAN provides precompiled binary versions of R for Windows, MacOS, and select Linux distributions that are likely sufficient for many users' needs.  Users can also install R from source code;  however this may require a significant amount of effort.  For specific details on how to compile, install, and manage R and R-packages, refer to the manual [R Installation and Administration](http://cran.r-project.org/doc/manuals/r-release/R-admin.html).

In its current construction, we recommend against running MEPIT while using R Studio.


### R Packages Required for MEPIT
MEPIT requires the installation of the following R libraries:

[doParallel](https://cran.r-project.org/web/packages/doParallel/index.html)

[Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html)

[RcppArmadillo](https://cran.r-project.org/web/packages/RcppArmadillo/index.html)

[RcppParallel](https://cran.r-project.org/web/packages/RcppParallel/index.html)

[CompQuadForm](https://cran.r-project.org/web/packages/CompQuadForm/index.html)

The easiest method to install these packages is with the following example command entered in an R shell:

    install.packages("doParallel", dependecies = TRUE)

Alternatively, one can also [install R packages from the command line]
              (http://cran.r-project.org/doc/manuals/r-release/R-admin.html#Installing-packages).

### Tutorial for Running MEPIT
For the simulation tutorial provided here, we generate genotypes for 3,000 samples typed at 10,000 unrelated variants. We show in our example R code how to implement MEPIT to perform a marginal epistasis association mapping test in order to find interacting causal variants of interest.

### Questions and Feedback
For questions or concerns with the MEPIT functions, please contact
[Lorin Crawford](mailto:lac55@stat.duke.edu) or 
[Xiang Zhou](mailto:xzhousph@umich.edu).

We appreciate any feedback you may have with our site and instructions.
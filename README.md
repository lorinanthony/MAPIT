# The MArginal ePIstasis Test (MAPIT)

### Introduction
Epistasis, commonly defined as the interaction between multiple genes, is an important genetic component underlying phenotypic variation. Many statistical methods have been developed to model and identify epistatic interactions between genetic variants. However, because of the large combinatorial search space of interactions, most epistasis mapping methods face enormous computational challenges and often suffer from low statistical power. In [Crawford et al. (2017)](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1006869), we present a novel, alternative strategy for mapping epistasis: the MArginal ePIstasis Test (MAPIT). Our method examines one variant at a time, and estimates and tests its "marginal epistatic effects" --- the combined pairwise interaction effects between a given variant and all other variants. By avoiding explicitly searching for interactions, our method avoids the large combinatorial search space and improves power. Our method is novel and relies on a recently developed variance component estimation method for efficient and robust parameter inference and p-value computation.

MAPIT is implemented as a set of R and C++ routines, which can be carried out within an R environment.


### The R Environment
R is a widely used, free, and open source software environment for statistical computing and graphics. The most recent version of R can be downloaded from the 
[Comprehensive R Archive Network (CRAN)](http://cran.r-project.org/)
CRAN provides precompiled binary versions of R for Windows, MacOS, and select Linux distributions that are likely sufficient for many users' needs.  Users can also install R from source code;  however, this may require a significant amount of effort.  For specific details on how to compile, install, and manage R and R-packages, refer to the manual [R Installation and Administration](http://cran.r-project.org/doc/manuals/r-release/R-admin.html).

In its current construction, we recommend against running MAPIT while using R Studio.


### R Packages Required for MAPIT
MAPIT requires the installation of the following R libraries:

[doParallel](https://cran.r-project.org/web/packages/doParallel/index.html)

[Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html)

[RcppArmadillo](https://cran.r-project.org/web/packages/RcppArmadillo/index.html)

[RcppParallel](https://cran.r-project.org/web/packages/RcppParallel/index.html)

[CompQuadForm](https://cran.r-project.org/web/packages/CompQuadForm/index.html)

The easiest method to install these packages is with the following example command entered in an R shell:

    install.packages("doParallel", dependecies = TRUE)

Alternatively, one can also [install R packages from the command line](http://cran.r-project.org/doc/manuals/r-release/R-admin.html#Installing-packages).

### C++ Functions Required for MAPIT
The code in this repository assumes that basic C++ functions and applications are already set up on the running personal computer or cluster. If not, the MAPIT functions and necessary Rcpp packages will not work properly. A simple option is to use [gcc](https://gcc.gnu.org/). macOS users may use this collection by installing the [Homebrew package manager](http://brew.sh/index.html) and then typing the following into the terminal:

    brew install gcc

For extra tips on how to run C++ on macOS, please visit [here](http://seananderson.ca/2013/11/18/rcpp-mavericks.html). For tips on how to avoid errors dealing with "-lgfortran" or "-lquadmath", please visit [here](http://thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/).

Note that there are two available versions of MAPIT, one of which takes advantage of [OpenMP](http://openmp.org/wp/), an API for multi-platform shared-memory parallel programming in C/C++. This is to speed up the computational time of the modeling algorithm. Unfortunately, OS X does not currently support OpenMP under the default compiler. A work around to use OpenMP in R on OS X can be found [here](http://thecoatlessprofessor.com/programming/openmp-in-r-on-os-x/). If after following these directions, the compiler of the C++ functions with 'omp.h' still does not work, we offer a standard version of MAPIT without the option of parallelization. If it is at all possible to get the OpenMP of MAPIT up and running, then we highly recommend it because this version of the code can drastically decrease computational time --- especially for larger data sets.
  

### Tutorial for Running MAPIT
For the simulation tutorial provided here, we generate genotypes for 3,000 samples typed at 10,000 unrelated variants. We show in our example R code how to implement MAPIT (both the standard and parallelized versions) to perform a marginal epistasis association mapping test in order to find interacting causal variants of interest.

### Relevant Citations
L. Crawford, P. Zeng, S. Mukherjee, X. Zhou (2017). Detecting epistasis with the marginal epistasis test in genetic mapping studies of quantitative traits. *PLoS Genet*. **13**(7): e1006869.

### Questions and Feedback
For questions or concerns with the MAPIT functions, please contact
[Lorin Crawford](mailto:lorin_crawford@brown.edu) or 
[Xiang Zhou](mailto:xzhousph@umich.edu).

We appreciate any feedback you may have with our repository and instructions.

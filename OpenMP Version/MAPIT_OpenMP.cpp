// load Rcpp
#include <RcppArmadillo.h>
#include <omp.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

// [[Rcpp::export]]
arma::mat GetLinearKernel(arma::mat X){
    double p = X.n_rows;
    return X.t()*X/p;
}

// [[Rcpp::export]]
arma::mat ComputePCs(arma::mat X,int top = 10){
    mat U;
    vec s;
    mat V;
    svd(U,s,V,X);
    
    mat PCs = U*diagmat(s);
    return PCs.cols(0,top-1);
}

////////////////////////////////////////////////////////////////////////////

//Below are functions for MAPIT using two hypothesis testing strategies:
//(1) MAPIT using the Normal or Z-Test
//(2) MAPIT using the Davies Method

//Considered are the following submodels:
//(1) Standard Model ---> y = m+g+e
//(2) Standard + Covariate Model ---> y = Wa+m+g+e
//(3) Standard + Common Environment Model ---> y = m+g+c+e
//(4) Standard + Covariate + Common Environment Model ---> y = Wa+m+g+c+e

//NOTE: delta = {delta(0),delta(1),delta(2)} = {sigma^2,omega^2,tau^2} for models (1) and (2)
//NOTE: delta = {delta(0),delta(1),delta(2),delta(3)} = {sigma^2,omega^2,nu^2,tau^2} for models (3) and (4)

////////////////////////////////////////////////////////////////////////////

// Model 1: Standard Model ---> y = m+g+e

// [[Rcpp::export]]
List MAPIT1_Normal(mat X,vec y,int cores = 1){
    int i;
    const int n = X.n_cols;
    const int p = X.n_rows;
    
    //Set up the vectors to save the outputs
    NumericVector sigma_est(p);
    NumericVector sigma_se(p);
    NumericVector pve(p);
    
    //Pre-compute the Linear GSM
    mat GSM = GetLinearKernel(X);
    
    omp_set_num_threads(cores);
#pragma omp parallel for schedule(dynamic)
    for(i=0; i<p; i++){
        
        //Compute K and G covariance matrices
        mat K = (GSM*p-trans(X.row(i))*X.row(i))/(p-1); //Create the linear kernel
        mat G = K; //Create the Kn Matrix
        G.each_row() %= X.row(i);
        G.each_col() %= trans(X.row(i));
        
        //Transform K and G using projection M
        mat b = zeros(n,2);
        b.col(0) = ones<vec>(n); b.col(1) = trans(X.row(i));
        mat btb_inv = inv(b.t()*b);
        mat Kc = K-b*btb_inv*(b.t()*K)-(K*b)*btb_inv*b.t()+b*btb_inv*(b.t()*(K*b))*btb_inv*b.t();
        mat Gc = G-b*btb_inv*(b.t()*G)-(G*b)*btb_inv*b.t()+b*btb_inv*(b.t()*(G*b))*btb_inv*b.t();
        vec yc = (eye<mat>(n,n)-(b*btb_inv)*b.t())*y;
        
        //Compute the quantities q and S
        vec q = zeros(3); //Create k-vector q to save
        mat S = zeros(3,3); //Create kxk-matrix S to save
        
        q(0) = as_scalar(yc.t()*Gc*yc);
        q(1) = as_scalar(yc.t()*Kc*yc);
        q(2) = as_scalar(yc.t()*(eye<mat>(n,n)-(b*btb_inv)*b.t())*yc);
        
        S(0,0) = as_scalar(accu(Gc.t()%Gc));
        S(0,1) = as_scalar(accu(Gc.t()%Kc));
        S(0,2) = as_scalar(accu(Gc.t()%(eye<mat>(n,n)-(b*btb_inv)*b.t())));
        S(1,0) = S(0,1);
        S(1,1) = as_scalar(accu(Kc.t()%Kc));
        S(1,2) = as_scalar(accu(Kc.t()%(eye<mat>(n,n)-(b*btb_inv)*b.t())));
        S(2,0) = S(0,2);
        S(2,1) = S(1,2);
        S(2,2) = as_scalar(accu(trans(eye<mat>(n,n)-(b*btb_inv)*b.t())%(eye<mat>(n,n)-(b*btb_inv)*b.t())));
        
        //Compute sigma(0) and sigma(1)
        mat Sinv = inv(S);
        vec delta = Sinv*q;
        
        //Compute var(sigma(0))
        double V_sigma = as_scalar(2*yc.t()*trans(Sinv(0,0)*Gc+Sinv(0,1)*Kc+Sinv(0,2)*(eye<mat>(n,n)-(b*btb_inv)*b.t()))*(delta(0)*Gc+delta(1)*Kc+delta(2)*(eye<mat>(n,n)-(b*btb_inv)*b.t()))*(Sinv(0,0)*Gc+Sinv(0,1)*Kc+Sinv(0,2)*(eye<mat>(n,n)-(b*btb_inv)*b.t()))*yc);
        
        //Save point estimates and SE of the epistasis component
        sigma_est(i) = delta(0);
        sigma_se(i) = sqrt(V_sigma);
        
        //Compute the PVE
        pve(i) = delta(0)/accu(delta);
    }
    
    //Compute the p-values for each estimate
    NumericVector sigma_pval = 2*Rcpp::pnorm(abs(sigma_est/sigma_se),0.0,1.0,0,0); //H0: sigma = 0 vs. H1: sigma != 0
    
    //Return a list of the arguments
    return Rcpp::List::create(Rcpp::Named("Est") = sigma_est, Rcpp::Named("SE") = sigma_se,Rcpp::Named("pvalues") = sigma_pval,Rcpp::Named("PVE") = pve);
}

////////////////////////////////////////////////////////////////////////////

// Model 2: Standard + Covariate Model ---> y = Wa+m+g+e

// [[Rcpp::export]]
List MAPIT2_Normal(mat X,vec y,mat Z,int cores = 1){
    int i;
    const int n = X.n_cols;
    const int p = X.n_rows;
    const int q = Z.n_rows;
    
    //Set up the vectors to save the outputs
    NumericVector sigma_est(p);
    NumericVector sigma_se(p);
    NumericVector pve(p);
    
    //Pre-compute the Linear GSM
    mat GSM = GetLinearKernel(X);
    
    omp_set_num_threads(cores);
#pragma omp parallel for schedule(dynamic)
    for(i=0; i<p; i++){
        
        //Compute K and G covariance matrices
        mat K = (GSM*p-trans(X.row(i))*X.row(i))/(p-1); //Create the linear kernel
        mat G = K; //Create the Kn Matrix
        G.each_row() %= X.row(i);
        G.each_col() %= trans(X.row(i));
        
        //Transform K and G using projection M
        mat b = zeros(n,q+2);
        b.col(0) = ones<vec>(n); b.cols(1,q) = Z.t(); b.col(q+1) = trans(X.row(i));
        mat btb_inv = inv(b.t()*b);
        mat Kc = K-b*btb_inv*(b.t()*K)-(K*b)*btb_inv*b.t()+b*btb_inv*(b.t()*(K*b))*btb_inv*b.t();
        mat Gc = G-b*btb_inv*(b.t()*G)-(G*b)*btb_inv*b.t()+b*btb_inv*(b.t()*(G*b))*btb_inv*b.t();
        vec yc = (eye<mat>(n,n)-(b*btb_inv)*b.t())*y;
        
        //Compute the quantities q and S
        vec q = zeros(3); //Create k-vector q to save
        mat S = zeros(3,3); //Create kxk-matrix S to save
        
        q(0) = as_scalar(yc.t()*Gc*yc);
        q(1) = as_scalar(yc.t()*Kc*yc);
        q(2) = as_scalar(yc.t()*(eye<mat>(n,n)-(b*btb_inv)*b.t())*yc);
        
        S(0,0) = as_scalar(accu(Gc.t()%Gc));
        S(0,1) = as_scalar(accu(Gc.t()%Kc));
        S(0,2) = as_scalar(accu(Gc.t()%(eye<mat>(n,n)-(b*btb_inv)*b.t())));
        S(1,0) = S(0,1);
        S(1,1) = as_scalar(accu(Kc.t()%Kc));
        S(1,2) = as_scalar(accu(Kc.t()%(eye<mat>(n,n)-(b*btb_inv)*b.t())));
        S(2,0) = S(0,2);
        S(2,1) = S(1,2);
        S(2,2) = as_scalar(accu(trans(eye<mat>(n,n)-(b*btb_inv)*b.t())%(eye<mat>(n,n)-(b*btb_inv)*b.t())));
        
        //Compute delta and Sinv
        mat Sinv = inv(S);
        vec delta = Sinv*q;
        
        //Compute var(delta(0))
        double V_sigma = as_scalar(2*yc.t()*trans(Sinv(0,0)*Gc+Sinv(0,1)*Kc+Sinv(0,2)*(eye<mat>(n,n)-(b*btb_inv)*b.t()))*(delta(0)*Gc+delta(1)*Kc+delta(2)*(eye<mat>(n,n)-(b*btb_inv)*b.t()))*(Sinv(0,0)*Gc+Sinv(0,1)*Kc+Sinv(0,2)*(eye<mat>(n,n)-(b*btb_inv)*b.t()))*yc);
        
        //Save point estimates and SE of the epistasis component
        sigma_est(i) = delta(0);
        sigma_se(i) = sqrt(V_sigma);
        
        //Compute the PVE
        pve(i) = delta(0)/accu(delta);
    }
    
    //Compute the p-values for each estimate
    NumericVector sigma_pval = 2*Rcpp::pnorm(abs(sigma_est/sigma_se),0.0,1.0,0,0); //H0: sigma = 0 vs. H1: sigma != 0
    
    //Return a list of the arguments
    return Rcpp::List::create(Rcpp::Named("Est") = sigma_est, Rcpp::Named("SE") = sigma_se,Rcpp::Named("pvalues") = sigma_pval,Rcpp::Named("PVE") = pve);
}

////////////////////////////////////////////////////////////////////////////

// Model 3: Standard + Common Environment Model ---> y = m+g+c+e

// [[Rcpp::export]]
List MAPIT3_Normal(mat X,vec y,mat C,int cores = 1){
    int i;
    const int n = X.n_cols;
    const int p = X.n_rows;
    
    //Set up the vectors to save the outputs
    NumericVector sigma_est(p);
    NumericVector sigma_se(p);
    NumericVector pve(p);
    
    //Pre-compute the Linear GSM
    mat GSM = GetLinearKernel(X);
    
    omp_set_num_threads(cores);
#pragma omp parallel for schedule(dynamic)
    for(i=0; i<p; i++){
        
        //Compute K and G covariance matrices
        mat K = (GSM*p-trans(X.row(i))*X.row(i))/(p-1); //Create the linear kernel
        mat G = K; //Create the Kn Matrix
        G.each_row() %= X.row(i);
        G.each_col() %= trans(X.row(i));
        
        //Transform K and G using projection M
        mat b = zeros(n,2);
        b.col(0) = ones<vec>(n); b.col(1) = trans(X.row(i));
        mat btb_inv = inv(b.t()*b);
        mat Kc = K-b*btb_inv*(b.t()*K)-(K*b)*btb_inv*b.t()+b*btb_inv*(b.t()*(K*b))*btb_inv*b.t();
        mat Gc = G-b*btb_inv*(b.t()*G)-(G*b)*btb_inv*b.t()+b*btb_inv*(b.t()*(G*b))*btb_inv*b.t();
        mat Cc = C-b*btb_inv*(b.t()*C)-(C*b)*btb_inv*b.t()+b*btb_inv*(b.t()*(C*b))*btb_inv*b.t();
        vec yc = (eye<mat>(n,n)-(b*btb_inv)*b.t())*y;
        
        //Compute the quantities q and S
        vec q = zeros(4); //Create k-vector q to save
        mat S = zeros(4,4); //Create kxk-matrix S to save
        
        q(0) = as_scalar(yc.t()*Gc*yc);
        q(1) = as_scalar(yc.t()*Kc*yc);
        q(2) = as_scalar(yc.t()*Cc*yc);
        q(3) = as_scalar(yc.t()*(eye<mat>(n,n)-(b*btb_inv)*b.t())*yc);
        
        S(0,0) = as_scalar(accu(Gc.t()%Gc));
        S(0,1) = as_scalar(accu(Gc.t()%Kc));
        S(0,2) = as_scalar(accu(Gc.t()%Cc));
        S(0,3) = as_scalar(accu(Gc.t()%(eye<mat>(n,n)-(b*btb_inv)*b.t())));
        
        S(1,0) = S(0,1);
        S(1,1) = as_scalar(accu(Kc.t()%Kc));
        S(1,2) = as_scalar(accu(Kc.t()%Cc));
        S(1,3) = as_scalar(accu(Kc.t()%(eye<mat>(n,n)-(b*btb_inv)*b.t())));
        
        S(2,0) = S(0,2);
        S(2,1) = S(1,2);
        S(2,2) = as_scalar(accu(Cc.t()%Cc));
        S(2,3) = as_scalar(accu(Cc.t()%(eye<mat>(n,n)-(b*btb_inv)*b.t())));
        
        S(3,0) = S(0,3);
        S(3,1) = S(1,3);
        S(3,2) = S(2,3);
        S(3,3) = as_scalar(accu(trans(eye<mat>(n,n)-(b*btb_inv)*b.t())%(eye<mat>(n,n)-(b*btb_inv)*b.t())));
        
        //Compute delta and Sinv
        mat Sinv = inv(S);
        vec delta = Sinv*q;
        
        //Compute var(delta(0))
        double V_sigma = as_scalar(2*yc.t()*trans(Sinv(0,0)*Gc+Sinv(0,1)*Kc+Sinv(0,2)*Cc+Sinv(0,3)*(eye<mat>(n,n)-(b*btb_inv)*b.t()))*(delta(0)*Gc+delta(1)*Kc+delta(2)*Cc+delta(3)*(eye<mat>(n,n)-(b*btb_inv)*b.t()))*(Sinv(0,0)*Gc+Sinv(0,1)*Kc+Sinv(0,2)*Cc+Sinv(0,3)*(eye<mat>(n,n)-(b*btb_inv)*b.t()))*yc);
        
        //Save point estimates and SE of the epistasis component
        sigma_est(i) = delta(0);
        sigma_se(i) = sqrt(V_sigma);
        
        //Compute the PVE
        pve(i) = delta(0)/accu(delta);
    }
    
    //Compute the p-values for each estimate
    NumericVector sigma_pval = 2*Rcpp::pnorm(abs(sigma_est/sigma_se),0.0,1.0,0,0); //H0: sigma = 0 vs. H1: sigma != 0
    
    //Return a list of the arguments
    return Rcpp::List::create(Rcpp::Named("Est") = sigma_est, Rcpp::Named("SE") = sigma_se,Rcpp::Named("pvalues") = sigma_pval,Rcpp::Named("PVE") = pve);
}

////////////////////////////////////////////////////////////////////////////

// Model 4: Standard + Covariate + Common Environment Model ---> y = Wa+m+g+c+e

// [[Rcpp::export]]
List MAPIT4_Normal(mat X,vec y,mat Z,mat C,int cores = 1){
    int i;
    const int n = X.n_cols;
    const int p = X.n_rows;
    const int q = Z.n_rows;
    
    //Set up the vectors to save the outputs
    NumericVector sigma_est(p);
    NumericVector sigma_se(p);
    NumericVector pve(p);
    
    //Pre-compute the Linear GSM
    mat GSM = GetLinearKernel(X);
    
    omp_set_num_threads(cores);
#pragma omp parallel for schedule(dynamic)
    for(i=0; i<p; i++){
        
        //Compute K and G covariance matrices
        mat K = (GSM*p-trans(X.row(i))*X.row(i))/(p-1); //Create the linear kernel
        mat G = K; //Create the Kn Matrix
        G.each_row() %= X.row(i);
        G.each_col() %= trans(X.row(i));
        
        //Transform K and G using projection M
        mat b = zeros(n,q+2);
        b.col(0) = ones<vec>(n); b.cols(1,q) = Z.t(); b.col(q+1) = trans(X.row(i));
        mat btb_inv = inv(b.t()*b);
        mat Kc = K-b*btb_inv*(b.t()*K)-(K*b)*btb_inv*b.t()+b*btb_inv*(b.t()*(K*b))*btb_inv*b.t();
        mat Gc = G-b*btb_inv*(b.t()*G)-(G*b)*btb_inv*b.t()+b*btb_inv*(b.t()*(G*b))*btb_inv*b.t();
        mat Cc = C-b*btb_inv*(b.t()*C)-(C*b)*btb_inv*b.t()+b*btb_inv*(b.t()*(C*b))*btb_inv*b.t();
        vec yc = (eye<mat>(n,n)-(b*btb_inv)*b.t())*y;
        
        //Compute the quantities q and S
        vec q = zeros(4); //Create k-vector q to save
        mat S = zeros(4,4); //Create kxk-matrix S to save
        
        q(0) = as_scalar(yc.t()*Gc*yc);
        q(1) = as_scalar(yc.t()*Kc*yc);
        q(2) = as_scalar(yc.t()*Cc*yc);
        q(3) = as_scalar(yc.t()*(eye<mat>(n,n)-(b*btb_inv)*b.t())*yc);
        
        S(0,0) = as_scalar(accu(Gc.t()%Gc));
        S(0,1) = as_scalar(accu(Gc.t()%Kc));
        S(0,2) = as_scalar(accu(Gc.t()%Cc));
        S(0,3) = as_scalar(accu(Gc.t()%(eye<mat>(n,n)-(b*btb_inv)*b.t())));
        
        S(1,0) = S(0,1);
        S(1,1) = as_scalar(accu(Kc.t()%Kc));
        S(1,2) = as_scalar(accu(Kc.t()%Cc));
        S(1,3) = as_scalar(accu(Kc.t()%(eye<mat>(n,n)-(b*btb_inv)*b.t())));
        
        S(2,0) = S(0,2);
        S(2,1) = S(1,2);
        S(2,2) = as_scalar(accu(Cc.t()%Cc));
        S(2,3) = as_scalar(accu(Cc.t()%(eye<mat>(n,n)-(b*btb_inv)*b.t())));
        
        S(3,0) = S(0,3);
        S(3,1) = S(1,3);
        S(3,2) = S(2,3);
        S(3,3) = as_scalar(accu(trans(eye<mat>(n,n)-(b*btb_inv)*b.t())%(eye<mat>(n,n)-(b*btb_inv)*b.t())));
        
        //Compute delta and Sinv
        mat Sinv = inv(S);
        vec delta = Sinv*q;
        
        //Compute var(delta(0))
        double V_sigma = as_scalar(2*yc.t()*trans(Sinv(0,0)*Gc+Sinv(0,1)*Kc+Sinv(0,2)*Cc+Sinv(0,3)*(eye<mat>(n,n)-(b*btb_inv)*b.t()))*(delta(0)*Gc+delta(1)*Kc+delta(2)*Cc+delta(3)*(eye<mat>(n,n)-(b*btb_inv)*b.t()))*(Sinv(0,0)*Gc+Sinv(0,1)*Kc+Sinv(0,2)*Cc+Sinv(0,3)*(eye<mat>(n,n)-(b*btb_inv)*b.t()))*yc);
        
        //Save point estimates and SE of the epistasis component
        sigma_est(i) = delta(0);
        sigma_se(i) = sqrt(V_sigma);
        
        //Compute the PVE
        pve(i) = delta(0)/accu(delta);
    }
    
    //Compute the p-values for each estimate
    NumericVector sigma_pval = 2*Rcpp::pnorm(abs(sigma_est/sigma_se),0.0,1.0,0,0); //H0: sigma = 0 vs. H1: sigma != 0
    
    //Return a list of the arguments
    return Rcpp::List::create(Rcpp::Named("Est") = sigma_est, Rcpp::Named("SE") = sigma_se,Rcpp::Named("pvalues") = sigma_pval,Rcpp::Named("PVE") = pve);
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

//Here is MAPIT using the Davies Method for hypothesis testing

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// Model 1: Standard Model ---> y = m+g+e

// [[Rcpp::export]]
List MAPIT1_Davies(mat X,vec y,const vec ind,int cores = 1){
    int i;
    const int n = X.n_cols;
    const int p = X.n_rows;
    
    //Set up the vectors to save the outputs
    NumericVector sigma_est(p);
    NumericVector pve(p);
    mat Lambda(n,p);
    
    //Pre-compute the Linear GSM
    mat GSM = GetLinearKernel(X);
    
    omp_set_num_threads(cores);
#pragma omp parallel for schedule(dynamic)
    for(i=0; i<p; i++){
        
        if(std::find(ind.begin(), ind.end(), i+1) != ind.end()){
            
            ///Compute K and G covariance matrices
            mat G = (GSM*p-trans(X.row(i))*X.row(i))/(p-1); //Create the linear kernel
            G.each_row() %= X.row(i);
            G.each_col() %= trans(X.row(i));
            
            //Transform K and G using projection M
            mat b = zeros(n,2);
            b.col(0) = ones<vec>(n); b.col(1) = trans(X.row(i));
            mat btb_inv = inv(b.t()*b);
            mat Gc = G-b*btb_inv*(b.t()*G)-(G*b)*btb_inv*b.t()+b*btb_inv*(b.t()*(G*b))*btb_inv*b.t(); //Scale Kn = MKnM^t
            vec yc = (eye<mat>(n,n)-(b*btb_inv)*b.t())*y;
            
            //Compute the quantities q and S
            vec q = zeros(2); //Create k-vector q to save
            mat S = zeros(2,2); //Create kxk-matrix S to save
            
            q(0) = as_scalar(yc.t()*Gc*yc);
            q(1) = as_scalar(yc.t()*(eye<mat>(n,n)-(b*btb_inv)*b.t())*yc);
            
            S(0,0) = as_scalar(accu(Gc.t()%Gc));
            S(0,1) = as_scalar(accu(Gc.t()%(eye<mat>(n,n)-(b*btb_inv)*b.t())));
            S(1,0) = S(0,1);
            S(1,1) = as_scalar(accu(trans(eye<mat>(n,n)-(b*btb_inv)*b.t())%(eye<mat>(n,n)-(b*btb_inv)*b.t())));
            
            //Compute delta and Sinv
            mat Sinv = inv(S);
            vec delta = Sinv*q;
            
            //Save point estimates and SE of the epistasis component
            sigma_est(i) = delta(0);
            
            //Find the eigenvalues of the projection matrix
            vec evals;
            eig_sym(evals,(Sinv(0,0)*Gc+Sinv(0,1)*(eye<mat>(n,n)-(b*btb_inv)*b.t()))*q(1)/S(1,1));
            Lambda.col(i) = evals;
            
            //Compute the PVE
            pve(i) = delta(0)/accu(delta);
        }
    }
    
    //Return a list of the arguments
    return Rcpp::List::create(Rcpp::Named("Est") = sigma_est, Rcpp::Named("Eigenvalues") = Lambda,Rcpp::Named("PVE") = pve);
}

////////////////////////////////////////////////////////////////////////////

// Model 2: Standard + Covariate Model ---> y = Wa+m+g+e

// [[Rcpp::export]]
List MAPIT2_Davies(mat X,vec y,mat Z,const vec ind,int cores = 1){
    int i;
    const int n = X.n_cols;
    const int p = X.n_rows;
    const int q = Z.n_rows;
    
    //Set up the vectors to save the outputs
    NumericVector sigma_est(p);
    NumericVector pve(p);
    mat Lambda(n,p);
    
    //Pre-compute the Linear GSM
    mat GSM = GetLinearKernel(X);
    
    omp_set_num_threads(cores);
#pragma omp parallel for schedule(dynamic)
    for(i=0; i<p; i++){
        
        if(std::find(ind.begin(), ind.end(), i+1) != ind.end()){
            
            //Compute K and G covariance matrices
            mat G = (GSM*p-trans(X.row(i))*X.row(i))/(p-1); //Create the linear kernel
            G.each_row() %= X.row(i);
            G.each_col() %= trans(X.row(i));
            
            //Transform K and G using projection M
            mat b = zeros(n,q+2);
            b.col(0) = ones<vec>(n); b.cols(1,q) = Z.t(); b.col(q+1) = trans(X.row(i));
            mat btb_inv = inv(b.t()*b);
            mat Gc = G-b*btb_inv*(b.t()*G)-(G*b)*btb_inv*b.t()+b*btb_inv*(b.t()*(G*b))*btb_inv*b.t();
            vec yc = (eye<mat>(n,n)-(b*btb_inv)*b.t())*y;
            
            //Compute the quantities q and S
            vec q = zeros(2); //Create k-vector q to save
            mat S = zeros(2,2); //Create kxk-matrix S to save
            
            q(0) = as_scalar(yc.t()*Gc*yc);
            q(1) = as_scalar(yc.t()*(eye<mat>(n,n)-(b*btb_inv)*b.t())*yc);
            
            S(0,0) = as_scalar(accu(Gc.t()%Gc));
            S(0,1) = as_scalar(accu(Gc.t()%(eye<mat>(n,n)-(b*btb_inv)*b.t())));
            S(1,0) = S(0,1);
            S(1,1) = as_scalar(accu(trans(eye<mat>(n,n)-(b*btb_inv)*b.t())%(eye<mat>(n,n)-(b*btb_inv)*b.t())));
            
            //Compute delta and Sinv
            mat Sinv = inv(S);
            vec delta = Sinv*q;
            
            //Save point estimates and SE of the epistasis component
            sigma_est(i) = delta(0);
            
            //Find the eigenvalues of the projection matrix
            vec evals;
            eig_sym(evals,(Sinv(0,0)*Gc+Sinv(0,1)*(eye<mat>(n,n)-(b*btb_inv)*b.t()))*q(1)/S(1,1));
            Lambda.col(i) = evals;
            
            //Compute the PVE
            pve(i) = delta(0)/accu(delta);
        }
    }
    
    //Return a list of the arguments
    return Rcpp::List::create(Rcpp::Named("Est") = sigma_est, Rcpp::Named("Eigenvalues") = Lambda,Rcpp::Named("PVE") = pve);
}

////////////////////////////////////////////////////////////////////////////

// Model 3: Standard + Common Environment Model ---> y = m+g+c+e

// [[Rcpp::export]]
List MAPIT3_Davies(mat X,vec y,mat C,const vec ind,int cores = 1){
    int i;
    const int n = X.n_cols;
    const int p = X.n_rows;
    
    //Set up the vectors to save the outputs
    NumericVector sigma_est(p);
    NumericVector pve(p);
    mat Lambda(n,p);
    
    //Pre-compute the Linear GSM
    mat GSM = GetLinearKernel(X);
    
    omp_set_num_threads(cores);
#pragma omp parallel for schedule(dynamic)
    for(i=0; i<p; i++){
        
        if(std::find(ind.begin(), ind.end(), i+1) != ind.end()){
            
            //Compute K and G covariance matrices
            mat K = (GSM*p-trans(X.row(i))*X.row(i))/(p-1); //Create the linear kernel
            mat G = K; //Create the Kn Matrix
            G.each_row() %= X.row(i);
            G.each_col() %= trans(X.row(i));
            
            //Transform K and G using projection M
            mat b = zeros(n,2);
            b.col(0) = ones<vec>(n); b.col(1) = trans(X.row(i));
            mat btb_inv = inv(b.t()*b);
            mat Kc = K-b*btb_inv*(b.t()*K)-(K*b)*btb_inv*b.t()+b*btb_inv*(b.t()*(K*b))*btb_inv*b.t();
            mat Gc = G-b*btb_inv*(b.t()*G)-(G*b)*btb_inv*b.t()+b*btb_inv*(b.t()*(G*b))*btb_inv*b.t();
            mat Cc = C-b*btb_inv*(b.t()*C)-(C*b)*btb_inv*b.t()+b*btb_inv*(b.t()*(C*b))*btb_inv*b.t();
            vec yc = (eye<mat>(n,n)-(b*btb_inv)*b.t())*y;
            
            //Compute the quantities q and S
            vec q = zeros(4); //Create k-vector q to save
            mat S = zeros(4,4); //Create kxk-matrix S to save
            
            q(0) = as_scalar(yc.t()*Gc*yc);
            q(1) = as_scalar(yc.t()*Kc*yc);
            q(2) = as_scalar(yc.t()*Cc*yc);
            q(3) = as_scalar(yc.t()*(eye<mat>(n,n)-(b*btb_inv)*b.t())*yc);
            
            S(0,0) = as_scalar(accu(Gc.t()%Gc));
            S(0,1) = as_scalar(accu(Gc.t()%Kc));
            S(0,2) = as_scalar(accu(Gc.t()%Cc));
            S(0,3) = as_scalar(accu(Gc.t()%(eye<mat>(n,n)-(b*btb_inv)*b.t())));
            
            S(1,0) = S(0,1);
            S(1,1) = as_scalar(accu(Kc.t()%Kc));
            S(1,2) = as_scalar(accu(Kc.t()%Cc));
            S(1,3) = as_scalar(accu(Kc.t()%(eye<mat>(n,n)-(b*btb_inv)*b.t())));
            
            S(2,0) = S(0,2);
            S(2,1) = S(1,2);
            S(2,2) = as_scalar(accu(Cc.t()%Cc));
            S(2,3) = as_scalar(accu(Cc.t()%(eye<mat>(n,n)-(b*btb_inv)*b.t())));
            
            S(3,0) = S(0,3);
            S(3,1) = S(1,3);
            S(3,2) = S(2,3);
            S(3,3) = as_scalar(accu(trans(eye<mat>(n,n)-(b*btb_inv)*b.t())%(eye<mat>(n,n)-(b*btb_inv)*b.t())));
            
            //Compute delta and Sinv
            mat Sinv = inv(S);
            vec delta = Sinv*q;
            
            //Record omega^2, nu^2, and tau^2 under the null hypothesis
            vec q_sub = zeros(3);
            mat S_sub = zeros(3,3);
            
            q_sub(0)=q(1);
            q_sub(1)=q(2);
            q_sub(2)=q(3);
            
            S_sub(0,0)=S(1,1);
            S_sub(0,1)=S(1,2);
            S_sub(0,2)=S(1,3);
            
            S_sub(1,0)=S(2,1);
            S_sub(1,1)=S(2,2);
            S_sub(1,2)=S(2,3);
            
            S_sub(2,0)=S(3,1);
            S_sub(2,1)=S(3,2);
            S_sub(2,2)=S(3,3);
            
            //Save point estimates and SE of the epistasis component
            sigma_est(i) = delta(0);
            
            //Compute P and P^{1/2} matrix
            vec delta_null = inv(S_sub)*q_sub;
            
            vec eigval;
            mat eigvec;
            
            eig_sym(eigval,eigvec,delta_null(0)*Kc+delta_null(1)*Cc+delta_null(2)*(eye<mat>(n,n)-(b*btb_inv)*b.t()));
            
            //Find the eigenvalues of the projection matrix
            vec evals;
            
            eig_sym(evals, (eigvec.cols(find(eigval>0))*diagmat(sqrt(eigval(find(eigval>0))))*trans(eigvec.cols(find(eigval>0))))*(Sinv(0,0)*Gc+Sinv(0,1)*Kc+Sinv(0,2)*Cc+Sinv(0,3)*(eye<mat>(n,n)-(b*btb_inv)*b.t()))*(eigvec.cols(find(eigval>0))*diagmat(sqrt(eigval(find(eigval>0))))*trans(eigvec.cols(find(eigval>0)))));
            Lambda.col(i) = evals;
            
            //Compute the PVE
            pve(i) = delta(0)/accu(delta);
        }
    }
    
    //Return a list of the arguments
    return Rcpp::List::create(Rcpp::Named("Est") = sigma_est, Rcpp::Named("Eigenvalues") = Lambda,Rcpp::Named("PVE") = pve);
}

////////////////////////////////////////////////////////////////////////////

// Model 4: Standard + Common Environment Model ---> y = Wa+m+g+c+e

// [[Rcpp::export]]
List MAPIT4_Davies(mat X,vec y,mat Z,mat C,const vec ind,int cores = 1){
    int i;
    const int n = X.n_cols;
    const int p = X.n_rows;
    const int q = Z.n_rows;
    
    //Set up the vectors to save the outputs
    NumericVector sigma_est(p);
    NumericVector pve(p);
    mat Lambda(n,p);
    
    //Pre-compute the Linear GSM
    mat GSM = GetLinearKernel(X);
    
    omp_set_num_threads(cores);
#pragma omp parallel for schedule(dynamic)
    for(i=0; i<p; i++){
        
        if(std::find(ind.begin(), ind.end(), i+1) != ind.end()){
            
            //Compute K and G covariance matrices
            mat K = (GSM*p-trans(X.row(i))*X.row(i))/(p-1); //Create the linear kernel
            mat G = K; //Create the Kn Matrix
            G.each_row() %= X.row(i);
            G.each_col() %= trans(X.row(i));
            
            //Transform K and G using projection M
            mat b = zeros(n,q+2);
            b.col(0) = ones<vec>(n); b.cols(1,q) = Z.t(); b.col(q+1) = trans(X.row(i));
            mat btb_inv = inv(b.t()*b);
            mat Kc = K-b*btb_inv*(b.t()*K)-(K*b)*btb_inv*b.t()+b*btb_inv*(b.t()*(K*b))*btb_inv*b.t();
            mat Gc = G-b*btb_inv*(b.t()*G)-(G*b)*btb_inv*b.t()+b*btb_inv*(b.t()*(G*b))*btb_inv*b.t();
            mat Cc = C-b*btb_inv*(b.t()*C)-(C*b)*btb_inv*b.t()+b*btb_inv*(b.t()*(C*b))*btb_inv*b.t();
            vec yc = (eye<mat>(n,n)-(b*btb_inv)*b.t())*y;
            
            //Compute the quantities q and S
            vec q = zeros(4); //Create k-vector q to save
            mat S = zeros(4,4); //Create kxk-matrix S to save
            
            q(0) = as_scalar(yc.t()*Gc*yc);
            q(1) = as_scalar(yc.t()*Kc*yc);
            q(2) = as_scalar(yc.t()*Cc*yc);
            q(3) = as_scalar(yc.t()*(eye<mat>(n,n)-(b*btb_inv)*b.t())*yc);
            
            S(0,0) = as_scalar(accu(Gc.t()%Gc));
            S(0,1) = as_scalar(accu(Gc.t()%Kc));
            S(0,2) = as_scalar(accu(Gc.t()%Cc));
            S(0,3) = as_scalar(accu(Gc.t()%(eye<mat>(n,n)-(b*btb_inv)*b.t())));
            
            S(1,0) = S(0,1);
            S(1,1) = as_scalar(accu(Kc.t()%Kc));
            S(1,2) = as_scalar(accu(Kc.t()%Cc));
            S(1,3) = as_scalar(accu(Kc.t()%(eye<mat>(n,n)-(b*btb_inv)*b.t())));
            
            S(2,0) = S(0,2);
            S(2,1) = S(1,2);
            S(2,2) = as_scalar(accu(Cc.t()%Cc));
            S(2,3) = as_scalar(accu(Cc.t()%(eye<mat>(n,n)-(b*btb_inv)*b.t())));
            
            S(3,0) = S(0,3);
            S(3,1) = S(1,3);
            S(3,2) = S(2,3);
            S(3,3) = as_scalar(accu(trans(eye<mat>(n,n)-(b*btb_inv)*b.t())%(eye<mat>(n,n)-(b*btb_inv)*b.t())));
            
            //Compute delta and Sinv
            mat Sinv = inv(S);
            vec delta = Sinv*q;
            
            //Record omega^2, nu^2, and tau^2 under the null hypothesis
            vec q_sub = zeros(3);
            mat S_sub = zeros(3,3);
            
            q_sub(0)=q(1);
            q_sub(1)=q(2);
            q_sub(2)=q(3);
            
            S_sub(0,0)=S(1,1);
            S_sub(0,1)=S(1,2);
            S_sub(0,2)=S(1,3);
            
            S_sub(1,0)=S(2,1);
            S_sub(1,1)=S(2,2);
            S_sub(1,2)=S(2,3);
            
            S_sub(2,0)=S(3,1);
            S_sub(2,1)=S(3,2);
            S_sub(2,2)=S(3,3);
            
            //Save point estimates and SE of the epistasis component
            sigma_est(i) = delta(0);
            
            //Compute P and P^{1/2} matrix
            vec delta_null = inv(S_sub)*q_sub;
            
            vec eigval;
            mat eigvec;
            
            eig_sym(eigval,eigvec,delta_null(0)*Kc+delta_null(1)*Cc+delta_null(2)*(eye<mat>(n,n)-(b*btb_inv)*b.t()));
            
            //Find the eigenvalues of the projection matrix
            vec evals;
            
            eig_sym(evals, (eigvec.cols(find(eigval>0))*diagmat(sqrt(eigval(find(eigval>0))))*trans(eigvec.cols(find(eigval>0))))*(Sinv(0,0)*Gc+Sinv(0,1)*Kc+Sinv(0,2)*Cc+Sinv(0,3)*(eye<mat>(n,n)-(b*btb_inv)*b.t()))*(eigvec.cols(find(eigval>0))*diagmat(sqrt(eigval(find(eigval>0))))*trans(eigvec.cols(find(eigval>0)))));
            Lambda.col(i) = evals;
            
            //Compute the PVE
            pve(i) = delta(0)/accu(delta);
        }
    }
    
    //Return a list of the arguments
    return Rcpp::List::create(Rcpp::Named("Est") = sigma_est, Rcpp::Named("Eigenvalues") = Lambda,Rcpp::Named("PVE") = pve);
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

//Here is MAPIT using the Davies Method for the hypothesis testing of every SNP

////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
List MAPIT_Davies(mat X,vec y,int cores = 1){
    int i;
    const int n = X.n_cols;
    const int p = X.n_rows;
    
    //Set up the vectors to save the outputs
    NumericVector sigma_est(p);
    NumericVector pve(p);
    mat Lambda(n,p);
    
    //Pre-compute the Linear GSM
    mat GSM = GetLinearKernel(X);
    
    omp_set_num_threads(cores);
#pragma omp parallel for schedule(dynamic)
    for(i=0; i<p; i++){
        ///Compute K and G covariance matrices
        mat G = (GSM*p-trans(X.row(i))*X.row(i))/(p-1); //Create the linear kernel
        G.each_row() %= X.row(i);
        G.each_col() %= trans(X.row(i));
        
        //Transform K and G using projection M
        mat b = zeros(n,2);
        b.col(0) = ones<vec>(n); b.col(1) = trans(X.row(i));
        mat btb_inv = inv(b.t()*b);
        mat Gc = G-b*btb_inv*(b.t()*G)-(G*b)*btb_inv*b.t()+b*btb_inv*(b.t()*(G*b))*btb_inv*b.t(); //Scale Kn = MKnM^t
        vec yc = (eye<mat>(n,n)-(b*btb_inv)*b.t())*y;
        
        //Compute the quantities q and S
        vec q = zeros(2); //Create k-vector q to save
        mat S = zeros(2,2); //Create kxk-matrix S to save
        
        q(0) = as_scalar(yc.t()*Gc*yc);
        q(1) = as_scalar(yc.t()*(eye<mat>(n,n)-(b*btb_inv)*b.t())*yc);
        
        S(0,0) = as_scalar(accu(Gc.t()%Gc));
        S(0,1) = as_scalar(accu(Gc.t()%(eye<mat>(n,n)-(b*btb_inv)*b.t())));
        S(1,0) = S(0,1);
        S(1,1) = as_scalar(accu(trans(eye<mat>(n,n)-(b*btb_inv)*b.t())%(eye<mat>(n,n)-(b*btb_inv)*b.t())));
        
        //Compute delta and Sinv
        mat Sinv = inv(S);
        vec delta = Sinv*q;
        
        //Save point estimates and SE of the epistasis component
        sigma_est(i) = delta(0);
        
        //Find the eigenvalues of the projection matrix
        vec evals;
        eig_sym(evals,(Sinv(0,0)*Gc+Sinv(0,1)*(eye<mat>(n,n)-(b*btb_inv)*b.t()))*q(1)/S(1,1));
        Lambda.col(i) = evals;
        
        //Compute the PVE
        pve(i) = delta(0)/accu(delta);
    }
    
    //Return a list of the arguments
    return Rcpp::List::create(Rcpp::Named("Est") = sigma_est, Rcpp::Named("Eigenvalues") = Lambda,Rcpp::Named("PVE") = pve);
}

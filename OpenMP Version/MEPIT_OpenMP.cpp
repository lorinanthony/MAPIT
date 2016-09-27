// load Rcpp
#include <RcppArmadillo.h>
#include <omp.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

// [[Rcpp::export]]
arma::mat GetLinearKernel(arma::mat X,int i){
    X.shed_row(i);
    double p = X.n_rows;
    return X.t()*X/p;
}

//Below are functions for MEPIT using two hypothesis testing strategies
//(1) MEPIT using the normal test
//(2) MEPIT using the Davies Method

//NOTE: {sigma(0),sigma(1),sigma(2)} = {sigma^2,omega^2,tau^2}

// [[Rcpp::export]]
List MEPIT_Normal(mat X,const vec y,int cores = 1){
    int i;
    const int n = X.n_cols;
    const int p = X.n_rows;
    
    //Set up the vectors to save the outputs
    NumericVector sigma_est(p);
    NumericVector sigma_se(p);
    NumericVector pve(p);
    
    omp_set_num_threads(cores);
#pragma omp parallel for schedule(dynamic)
    for(i=0; i<p; i++){
        
        //Compute K and G covariance matrices
        mat K = GetLinearKernel(X,i); //Create the linear kernel
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
        vec sigma = zeros(3);
        mat Sinv = inv(S);
        sigma = Sinv*q;
        
        //Compute var(sigma(0))
        double V_sigma = as_scalar(2*yc.t()*trans(Sinv(0,0)*Gc+Sinv(0,1)*Kc+Sinv(0,2)*(eye<mat>(n,n)-(b*btb_inv)*b.t()))*(sigma(0)*Gc+sigma(1)*Kc+sigma(2)*(eye<mat>(n,n)-(b*btb_inv)*b.t()))*(Sinv(0,0)*Gc+Sinv(0,1)*Kc+Sinv(0,2)*(eye<mat>(n,n)-(b*btb_inv)*b.t()))*yc);
        
        //Save point estimates and SE of the epistasis component
        sigma_est(i) = sigma(0);
        sigma_se(i) = sqrt(V_sigma);
        
        //Compute the PVE
        pve(i) = sigma(0)/accu(sigma);
    }
    
    //Compute the p-values for each estimate
    NumericVector sigma_pval = 2*Rcpp::pnorm(abs(sigma_est/sigma_se),0.0,1.0,0,0); //H0: sigma = 0 vs. H1: sigma != 0
    
    //Return a list of the arguments
    return Rcpp::List::create(Rcpp::Named("Est") = sigma_est, Rcpp::Named("SE") = sigma_se,Rcpp::Named("pvalues") = sigma_pval,Rcpp::Named("PVE") = pve);

}

// [[Rcpp::export]]
List MEPIT_Davies(mat X,const vec y,const vec ind,int cores = 1){
    int i;
    const int n = X.n_cols;
    const int p = X.n_rows;
    
    //Set up the vectors to save the outputs
    NumericVector sigma_est(p);
    NumericVector pve(p);
    mat Lambda(n,p);
    
    omp_set_num_threads(cores);
#pragma omp parallel for schedule(dynamic)
    for(i=0; i<p; i++){
        
        if(std::find(ind.begin(), ind.end(), i+1) != ind.end()){
            
            //Compute K and G covariance matrices
            mat K = GetLinearKernel(X,i); //Create the linear kernel
            mat G = K; //Create the Kn Matrix
            G.each_row() %= X.row(i);
            G.each_col() %= trans(X.row(i));
            
            //Transform K and G using projection M
            mat b = zeros(n,2);
            b.col(0) = ones<vec>(n); b.col(1) = trans(X.row(i));
            mat btb_inv = inv(b.t()*b);
            mat Kc = K-b*btb_inv*(b.t()*K)-(K*b)*btb_inv*b.t()+b*btb_inv*(b.t()*(K*b))*btb_inv*b.t();
            mat Gc = G-b*btb_inv*(b.t()*G)-(G*b)*btb_inv*b.t()+b*btb_inv*(b.t()*(G*b))*btb_inv*b.t(); //Scale Kn = MKnM^t
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
            
            //Record omega^2 and tau^2 under the null hypothesis
            vec q_sub = zeros(2);
            mat S_sub = zeros(2,2);
            
            q_sub(0)=q(1);
            q_sub(1)=q(2);
            
            S_sub(0,0)=S(1,1);
            S_sub(0,1)=S(1,2);
            S_sub(1,0)=S(2,1);
            S_sub(1,1)=S(2,2);
            
            //Compute sigma(0) and sigma(1)
            vec sigma = zeros(3);
            mat Sinv = inv(S);
            sigma = Sinv*q;
            
            //Save point estimates and SE of the epistasis component
            sigma_est(i) = sigma(0);
            
            //Compute P and P^{1/2} matrix
            vec sigma_null = zeros(2);
            sigma_null = inv(S_sub)*q_sub;
            
            vec eigval;
            mat eigvec;
            
            eig_sym(eigval,eigvec,sigma_null(0)*Kc+sigma_null(1)*(eye<mat>(n,n)-(b*btb_inv)*b.t()));
            
            //Find the eigenvalues of the projection matrix
            vec evals;
            
            eig_sym(evals, (eigvec.cols(find(eigval>0))*diagmat(sqrt(eigval(find(eigval>0))))*trans(eigvec.cols(find(eigval>0))))*(Sinv(0,0)*Gc+Sinv(0,1)*Kc+Sinv(0,2)*(eye<mat>(n,n)-(b*btb_inv)*b.t()))*(eigvec.cols(find(eigval>0))*diagmat(sqrt(eigval(find(eigval>0))))*trans(eigvec.cols(find(eigval>0)))));
            Lambda.col(i) = evals;
            
            //Compute the PVE
            pve(i) = sigma(0)/accu(sigma);
        }
    }
    
    //Return a list of the arguments
    return Rcpp::List::create(Rcpp::Named("Est") = sigma_est, Rcpp::Named("Eigenvalues") = Lambda,Rcpp::Named("PVE") = pve);
}
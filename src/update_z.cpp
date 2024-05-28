// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector update_z(NumericVector zs, 
                       List A, 
                       NumericMatrix Ps,
                       NumericVector pis,
                       NumericVector classes) 
{
    int L = A.length();
    int n = zs.length();
    int K0 = Ps.ncol();
    NumericVector z_ret(n);
    for(int i = 0; i < n; i ++)
    {
        NumericVector pi_star(K0);
        for(int k = 0; k < K0; k++)
        {
            pi_star[k] = log(pis[k]);
            for(int j = 0; j < n; j++)
            {
                if(i != j)
                {
                    for(int l = 0; l < L; l++)
                    {
                        NumericMatrix Al = A[l];
                        pi_star[k] = ((pi_star[k]) + 
                            log(pow(Ps(k,zs[j]-1),Al(i,j))) + 
                            log(pow((1 - Ps(k,zs[j]-1)),(1-Al(i,j)))));
                    }
                }
            }
            //Rcout << i << "," << k << ":" << pi_star[k] << std::endl;
        }
        //Rcout << i << ":" << pi_star << std::endl;
        NumericVector pi_star2 (pi_star.length());
        if(any(abs(pi_star) > 700).is_true())
        {
            pi_star2[which_max(pi_star)] = 1;
        }
        else
        {
            pi_star2 = exp(pi_star) / sum(exp(pi_star));
        }
        
        // pi_star2 = exp(pi_star) / sum(exp(pi_star));
        //Rcout << i << ":" << pi_star2 << std::endl;
        NumericVector zi = Rcpp::RcppArmadillo::sample(classes,1,TRUE,pi_star2);
        z_ret[i] = zi[0];
    }
    return z_ret;
}


// [[Rcpp::export]]
NumericVector update_z_single(NumericVector zs, NumericMatrix A, NumericMatrix Ps,NumericVector pis,NumericVector classes) 
{
    int n = zs.length();
    int K0 = Ps.ncol();
    NumericVector z_ret(n);
    for(int i = 0; i < n; i ++)
    {
        NumericVector pi_star(K0);
        for(int k = 0; k < K0; k++)
        {
            pi_star[k] = log(pis[k]);
            for(int j = 0; j < n; j++)
            {
                if(i != j)
                {
                    // pi_star[k] = pi_star[k] * 
                    //     pow(Ps(k,zs[j]-1),A(i,j)) * 
                    //     pow((1 - Ps(k,zs[j]-1)),(1-A(i,j)));
                    pi_star[k] = ((pi_star[k]) + 
                        log(pow(Ps(k,zs[j]-1),A(i,j))) + 
                        log(pow((1 - Ps(k,zs[j]-1)),(1-A(i,j)))));
                }
            }
            //Rcout << i << "," << k << ":" << pi_star[k] << std::endl;
        }
        //pi_star = pi_star / n;
        //Rcout << i << ":" << pi_star << std::endl;
        NumericVector pi_star2 (pi_star.length());
        if(any(abs(pi_star) > 700).is_true())
        {
            pi_star2[which_max(pi_star)] = 1;
        }
        else
        {
            pi_star2 = exp(pi_star) / sum(exp(pi_star));
        }

        // pi_star2 = exp(pi_star) / sum(exp(pi_star));
        //Rcout << i << ":" << pi_star2 << std::endl;
        NumericVector zi = Rcpp::RcppArmadillo::sample(classes,1,TRUE,pi_star2);
        z_ret[i] = zi[0];
    }
    return z_ret;
}

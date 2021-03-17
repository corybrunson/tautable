#include <Rcpp.h>
include namespace Rcpp;

// [[Rcpp::export]]
NumericVector mahonians(int n) {
    NumericVector vec = 1;
    
    for(int m = 1; m < n; m++) {
        int p = vec.size();
        NumericMatrix mat(m + 1, p + m);
        for(int i = 0; i < m + 1, i++) {
            for(int j = 0; j < p + m; j++) {
                
            }
        }
        
        NumericVector vec(nrow)
        
        for(int i = 0; i < m + 1; i++) {
            double total = 0;
            for(int j = 0; j < p + m; j++) {
                total += mat(i, j);
            }
            vec[i] = total;
        }
    }
}

// [[Rcpp::export]]
NumericVector mahonians(int n) {
    
    int l = 1 + n * (n - 1) / 2;
    NumericVector x(l);
    x[1] = 1;
    
    for (int m = 1; m <= n; m++) {
        
        int p = x.size();
        NumericVector y(p + m);
        
        
        
        x = y
    }
}
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix psiFastCalc(const NumericMatrix& mat,
                          const NumericVector incA=0, 
                          const NumericVector incB=0,
                          const NumericVector excA=0, 
                          const NumericVector excB=0,
                          const int minReads=10) {
    double inc, exc, total, psi;
    NumericMatrix out(incA.length(), mat.ncol());
    
    for (size_t col=0; col < mat.ncol(); col++) {
        for (size_t idx=0; idx < incA.length(); idx++) {
            inc = mat(incA[idx] - 1, col);
            if ( incB[0] > 0 ) inc = (inc + mat(incB[idx] - 1, col))/2;
            
            exc = mat(excA[idx] - 1, col);
            if ( excB[0] > 0 ) exc = (exc + mat(excB[idx] - 1, col))/2;
            
            total = inc + exc;
            if (total < minReads)
                psi = NumericVector::get_na();
            else
                psi = inc / total;
            out(idx, col) = psi;
        }
    }
    colnames(out) = colnames(mat);
    return out;
}

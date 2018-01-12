#include <Rcpp.h>
using namespace Rcpp;

void progressBar(double progress) {
    // Source: http://lists.r-forge.r-project.org/pipermail/rcpp-devel/2010-August/000964.html
    
    // Create progress bar
    int barWidth=40;
    printf("  |");
    
    // Print completed progress
    int complete = round(progress * barWidth);
    for (int i=0; i < complete; i++) printf("=");
    // Print remaining progress
    for (int i=complete; i < barWidth; i++) printf(" ");
    printf("| %3.0f%% \r", progress * 100);
    fflush(stdout); // Avoids output buffering problems
}

// [[Rcpp::export]]
NumericMatrix psiFastCalc(const NumericMatrix& mat,
                          const NumericVector incA=0, 
                          const NumericVector incB=0,
                          const NumericVector excA=0, 
                          const NumericVector excB=0,
                          const int minReads=10) {
    double inc, exc, total, psi, progress;
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
        progress = col / (mat.ncol() - 1);
        progressBar(progress);
    }
    colnames(out) = colnames(mat);
    return out;
}

// [[Rcpp::export]]
NumericMatrix psiFastCalc2(const NumericMatrix& mat,
                           const List& inc, const List& exc, 
                           const int minReads=10) {
    double incReads, excReads, totalReads, psi, progress;
    NumericMatrix out(inc.length(), mat.ncol());
    
    for (size_t col=0; col < mat.ncol(); col++) {
        for (size_t idx=0; idx < inc.length(); idx++) {
            NumericVector incIdx = as<NumericVector>(inc[idx]);
            incReads = 0;
            for (size_t k=0; k < incIdx.length(); k++)
                incReads = incReads + mat(incIdx[k] - 1, col);
            
            NumericVector excIdx = as<NumericVector>(exc[idx]);
            excReads = 0;
            for (size_t k=0; k < excIdx.length(); k++)
                excReads = excReads + mat(excIdx[k] - 1, col);
            
            totalReads = incReads + excReads;
            if (totalReads < minReads)
                psi = NumericVector::get_na();
            else
                psi = incReads / totalReads;
            out(idx, col) = psi;
        }
        progress = col / (mat.ncol() - 1);
        progressBar(progress);
    }
    colnames(out) = colnames(mat);
    return out;
}

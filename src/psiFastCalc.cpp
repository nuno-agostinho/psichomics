#include <Rcpp.h>
using namespace Rcpp;

void progressBar (double progress) {
    // Source: https://lists.r-forge.r-project.org/pipermail/rcpp-devel/2010-August/000964.html
    
    // Create progress bar
    int barWidth=40;
    Rprintf("  |");
    
    // Print completed progress
    int complete = round(progress * barWidth);
    for (int i = 0; i < complete; i++) Rprintf("=");
    // Print remaining progress
    for (int i = complete; i < barWidth; i++) Rprintf(" ");
    Rprintf("| %3.0f%% \r", progress * 100);
    if (progress == 1) Rprintf("\n");
}

double calculatePSI (double inc, double exc, double minReads) {
    double psi, total = inc + exc;
    if (total < minReads) {
        psi = NumericVector::get_na();
    } else {
        psi = inc / total;
    }
    return psi;
}

// [[Rcpp::export]]
NumericMatrix psiFastCalc(const NumericMatrix& mat,
                          const NumericVector incA=0, 
                          const NumericVector incB=0,
                          const NumericVector excA=0, 
                          const NumericVector excB=0,
                          const int minReads=10) {
    double incReads, excReads;
    int ncol = mat.ncol(), incLen = incA.length();
    NumericMatrix out(incLen, ncol);
    
    for (int col = 0; col < ncol; col++) {
        for (int idx = 0; idx < incLen; idx++) {
            incReads = mat(incA[idx] - 1, col);
            if ( incB[0] > 0 ) {
                incReads = (incReads + mat(incB[idx] - 1, col))/2;
            }
            excReads = mat(excA[idx] - 1, col);
            if ( excB[0] > 0 ) {
                excReads = (excReads + mat(excB[idx] - 1, col))/2;
            }
            out(idx, col) = calculatePSI(incReads, excReads, minReads);
        }
        if (ncol > 1) progressBar(double(col) / double(mat.ncol() - 1));
    }
    colnames(out) = colnames(mat);
    return out;
}

// [[Rcpp::export]]
NumericMatrix psiFastCalc2(const NumericMatrix& mat,
                           const List& inc, const List& exc, 
                           const int minReads=10) {
    double incReads, excReads;
    int ncol = mat.ncol(), incLen = inc.length();
    NumericMatrix out(incLen, ncol);
    NumericVector incIdx, excIdx;
    
    for (int col = 0; col < ncol; col++) {
        for (int idx = 0; idx < incLen; idx++) {
            incIdx = as<NumericVector>(inc[idx]);
            excIdx = as<NumericVector>(exc[idx]);
            incReads = 0;
            for (int k = 0; k < incIdx.length(); k++) {
                incReads = incReads + mat(incIdx[k] - 1, col);
            }
            excReads = 0;
            for (int k = 0; k < excIdx.length(); k++) {
                excReads = excReads + mat(excIdx[k] - 1, col);
            }
            out(idx, col) = calculatePSI(incReads, excReads, minReads);
        }
        if (mat.ncol() > 1) progressBar(double(col) / double(mat.ncol() - 1));
    }
    colnames(out) = colnames(mat);
    return out;
}

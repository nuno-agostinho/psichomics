#include <Rcpp.h>
#include "progressBar.h"
using namespace Rcpp;

bool isIndicatedCvg (CharacterVector cvg, int row,
                     CharacterVector scoresToDiscard) {
    // Check if coverage/quality values are the ones to discard
    auto quality = cvg[row];
    bool isToDiscard = false;
    for (auto score : scoresToDiscard) {
        int l = 0;
        for (int k = 0; k < quality.size(); k++) {
            if (quality[k] == score[l]) {
                // While characters are matching
                l++;
                if (l >= score.size()) {
                    isToDiscard = true;
                    break;
                }
            } else if (l > 0) {
                // If first character but not subsequent ones are found
                break;
            }
        }
        if (isToDiscard) break;
    }
    return isToDiscard;
}

// [[Rcpp::export]]
DataFrame discardVastToolsByCvg(DataFrame psi, DataFrame eventData,
                                int qualityCol,
                                CharacterVector scoresToDiscard) {
    DataFrame res = Rcpp::clone(psi);
    // Prepare scores to be matched unequivocally
    CharacterVector scores = Rcpp::clone(scoresToDiscard);
    for (int n = 0; n < scores.length(); n++) scores[n] = "," + scores[n] + ",";

    int ncol = psi.length(), nrow = psi.nrows();
    CharacterVector sampleCvg;
    NumericVector sampleVal;
    for (int col = 0; col < ncol; col++) {
        sampleVal = res[col];
        sampleCvg = eventData[col + qualityCol];
        for (int row = 0; row < nrow; row++) {
            if ( isIndicatedCvg(sampleCvg, row, scores) ) {
                sampleVal[row] = NA_REAL;
            }
        }
        progressBar(double(col) / double(ncol - 1));
    }
    return res;
}

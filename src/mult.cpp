#include <Rcpp.h>
using namespace Rcpp;

//' @export
// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  return x * 2;
}

// Simple testing
/*** R
timesTwo(42)
*/

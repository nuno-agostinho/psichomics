#include <Rcpp.h>
using namespace Rcpp;

void progressBar (double progress) {
    // Source:
    // https://lists.r-forge.r-project.org/pipermail/rcpp-devel/2010-August/000964.html

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

#include "deoverlap.h"
#include <algorithm>

using namespace Rcpp;

RcppExport SEXP deoverlap(SEXP start_, SEXP end_) {
    NumericVector start(start_);
    NumericVector end(end_);

    for (int i = 0; i < start.size() - 1; ++i) {
        int j = i + 1;
        if (end[i] > start[j]) {
            double midpoint = (end[i] + start[j]) / 2;
            end[i] = midpoint;
            start[j] = midpoint;
        }
    }

    return List::create(
        Named("start") = start,
        Named("end") = end);
}

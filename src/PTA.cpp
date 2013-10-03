#include "PTA.h"
#include "PTAProcessor.h"
using namespace Rcpp;
using namespace std;

SEXP PTA(SEXP start, SEXP end, SEXP score,
        SEXP count_bound, SEXP error_bound, SEXP adjacency_threshold, SEXP skip, SEXP mode,
        SEXP correlation_bound, SEXP correlation_spearman) {
BEGIN_RCPP
    PTAProcessor p(start, end, score, count_bound, error_bound, adjacency_threshold, skip, mode, correlation_bound, correlation_spearman);
    return p.run();
END_RCPP
}

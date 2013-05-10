#include "PTA.h"
#include "PTAProcessor.h"
using namespace Rcpp;
using namespace std;

SEXP PTA(SEXP start, SEXP end, SEXP score, SEXP count, SEXP error_bound, SEXP adjacency_treshold, SEXP skip, SEXP mode, SEXP correlation_bound) {
BEGIN_RCPP
    PTAProcessor p(start, end, score, count, error_bound, adjacency_treshold, skip, mode, correlation_bound);
    return p.run();
END_RCPP
}

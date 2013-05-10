#include "PTA.h"
#include "PTAProcessor.h"
using namespace Rcpp;
using namespace std;

SEXP PTA(SEXP start, SEXP end, SEXP score, SEXP count, SEXP error, SEXP adjacency_treshold, SEXP skip, SEXP mode, SEXP correlation_bound) {
BEGIN_RCPP
    PTAProcessor p(start, end, score, count, error, adjacency_treshold, skip, mode, correlation_bound);
    p.run();
    return p.get_result();
END_RCPP
}

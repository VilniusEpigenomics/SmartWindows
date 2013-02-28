#include "pta.h"
#include "OptimalPTAProcessor.h"
#include "GreedyPTAProcessor.h"
using namespace Rcpp;
using namespace std;

SEXP PTA(SEXP start, SEXP end, SEXP score, SEXP count, SEXP error, SEXP adjacency_treshold) {
BEGIN_RCPP
    OptimalPTAProcessor p(start, end, score, count, error, adjacency_treshold);
    p.run();
    return p.get_result();
END_RCPP
}

SEXP gPTA(SEXP start, SEXP end, SEXP score, SEXP count, SEXP error, SEXP adjacency_treshold) {
BEGIN_RCPP
    GreedyPTAProcessor p(start, end, score, count, error, adjacency_treshold);
    p.run();
    return p.get_result();
END_RCPP
}

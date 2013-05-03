#include "pta.h"
#include "OptimalPTAProcessor.h"
#include "GreedyPTAProcessor.h"
#include "MultiPTAProcessor.h"
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

SEXP multiPTA(SEXP start, SEXP end, SEXP score, SEXP count, SEXP error, SEXP adjacency_treshold, SEXP skip) {
BEGIN_RCPP
    MultiPTAProcessor p(start, end, score, count, error, adjacency_treshold, skip);
    p.run();
    return p.get_result();
END_RCPP
}

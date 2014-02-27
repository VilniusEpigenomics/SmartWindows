#include "aggregators.h"
#include "PTAProcessor.h"
#include "SpanAggregator.h"
using namespace Rcpp;
using namespace std;

SEXP PTA(SEXP arguments_) {
BEGIN_RCPP
    const List arguments(arguments_);
    PTAProcessor p(arguments);
    return p.run();
END_RCPP
}

SEXP spanAggregate(SEXP arguments_) {
BEGIN_RCPP
    const List arguments(arguments_);
    SpanAggregator a(arguments);
    return a.run();
END_RCPP
}

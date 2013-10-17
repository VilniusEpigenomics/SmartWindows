#include "PTA.h"
#include "PTAProcessor.h"
using namespace Rcpp;
using namespace std;

SEXP PTA(SEXP arguments_) {
BEGIN_RCPP
    const List arguments(arguments_);
    PTAProcessor p(arguments);
    return p.run();
END_RCPP
}

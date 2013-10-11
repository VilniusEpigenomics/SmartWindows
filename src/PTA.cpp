#include "PTA.h"
#include "PTAProcessor.h"
using namespace Rcpp;
using namespace std;

SEXP PTA(SEXP arguments) {
BEGIN_RCPP
    PTAProcessor p(arguments);
    return p.run();
END_RCPP
}

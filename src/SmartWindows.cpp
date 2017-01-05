#include "SmartWindowsProcessor.h"
using namespace Rcpp;

RcppExport SEXP SmartWindows(SEXP arguments_) {
BEGIN_RCPP
    const List arguments(arguments_);
    SmartWindowsProcessor p(arguments);
    return p.run();
END_RCPP
}

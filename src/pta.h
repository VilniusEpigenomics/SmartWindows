#ifndef pta_pta_h
#define pta_pta_h

#include <Rcpp.h>

RcppExport __attribute__((visibility("default"))) SEXP PTA(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport __attribute__((visibility("default"))) SEXP gPTA(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport __attribute__((visibility("default"))) SEXP multiPTA(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

#endif

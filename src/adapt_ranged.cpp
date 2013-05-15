#include <Rcpp.h>
#include "adapt_ranged.h"

using namespace Rcpp;

RcppExport SEXP adapt_ranged(
        SEXP sum_, SEXP len_, SEXP found_,
        SEXP error_, SEXP add_error_,
        SEXP destnrow_, SEXP deststart_, SEXP destend_, SEXP destscore_,
        SEXP srcnrow_, SEXP srcstart_, SEXP srcend_, SEXP srcscore_)
{
    NumericVector sum(sum_);
    NumericVector len(len_);
    LogicalVector found(found_);
    NumericVector error(error_);
    int add_error = as<int>(add_error_);
    int destnrow = as<int>(destnrow_);
    NumericVector deststart(deststart_);
    NumericVector destend(destend_);
    NumericVector destscore(destscore_);
    int srcnrow = as<int>(srcnrow_);
    NumericVector srcstart(srcstart_);
    NumericVector srcend(srcend_);
    NumericVector srcscore(srcscore_);

    int srci = 0;
    int desti = 0;
    double intlen, left, right;
    bool rightinside; // src inside dest

    while ((srci < srcnrow) && (desti < destnrow)) {
        if (deststart[desti] >= srcend[srci]) {
            ++srci;
            continue;
        }

        if (srcstart[srci] >= destend[desti]) {
            ++desti;
            continue;
        }

        if (srcstart[srci] >= deststart[desti]) {
            left = srcstart[srci];
        } else {
            left = deststart[desti];
        }

        if (srcend[srci] <= destend[desti]) {
            right = srcend[srci];
            rightinside = true;
        } else {
            right = destend[desti];
            rightinside = false;
        }

        intlen = right - left;
        found[desti] = true;
        sum[desti] += intlen * srcscore[srci];
        len[desti] += intlen;
        if (add_error) error[desti] += intlen * pow(srcscore[srci] - destscore[desti], 2);

        if (rightinside) {
            ++srci;
        } else {
            ++desti;
        }
    }

    return R_NilValue;
}

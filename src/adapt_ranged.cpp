#include <Rcpp.h>
#include "adapt_ranged.h"

using namespace Rcpp;

RcppExport SEXP adapt_ranged(
    SEXP srcstart_, SEXP srcend_, SEXP srcscore_,
    SEXP deststart_, SEXP destend_,
    SEXP add_error_)
{
    NumericVector srcstart(srcstart_);
    NumericVector srcend(srcend_);
    NumericVector srcscore(srcscore_);
    NumericVector deststart(deststart_);
    NumericVector destend(destend_);
    int add_error = as<int>(add_error_);

    int srcn = srcstart.size();
    int destn = deststart.size();

    NumericVector len(destn, 0.0);
    NumericVector sum(destn, 0.0);
    LogicalVector found(destn, false);
    NumericVector error(destn, 0.0);
    NumericVector destscore(destn, 0.0);

    int srci = 0;
    int desti = 0;
    double intlen, left, right;
    bool rightinside; // src inside dest

    while ((srci < srcn) && (desti < destn)) {
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

    return wrap(ifelse(found, sum / len, NA_REAL));
}

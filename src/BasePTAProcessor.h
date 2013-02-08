#ifndef pta_BasePTAProcessor_h
#define pta_BasePTAProcessor_h
#include <Rcpp.h>

class BasePTAProcessor {
    protected:
        Rcpp::NumericVector start;
        Rcpp::NumericVector end;
        Rcpp::NumericVector score;

        int count_bound;
        double error_bound;
        bool error_bounded;

        int minimum_count;
        double maximum_error;

        inline int size() const { return start.size(); }
        double length(int interval) const;
        bool adjacent(int i, int j) const;
        double merged_score(int i, int j) const;

    public:
        BasePTAProcessor(SEXP start_, SEXP end_, SEXP score_, SEXP count_, SEXP error_);
        Rcpp::List get_result() const;
};

#endif

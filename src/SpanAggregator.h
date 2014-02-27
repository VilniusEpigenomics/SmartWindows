#ifndef PTA_SpanAggregator_h
#define PTA_SpanAggregator_h
#include <Rcpp.h>
#include <vector>
#include <cassert>

class __attribute__((visibility("hidden"))) SpanAggregator {
    private:
        const Rcpp::IntegerVector start;
        const Rcpp::IntegerVector end;
        const Rcpp::NumericMatrix scores;

        const long span;

        inline int size() const { return start.size(); }

    public:
        SpanAggregator(const Rcpp::List arguments);
        Rcpp::List run();
};

#endif

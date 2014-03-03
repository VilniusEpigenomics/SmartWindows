#ifndef PTA_IntersectionAggregator_h
#define PTA_IntersectionAggregator_h
#include <Rcpp.h>
#include <deque>
#include <cassert>

class IntersectionAggregator {
    private:
        struct Range {
            long group;
            long start;
            long end;
            Rcpp::NumericVector score;
        };

        const Rcpp::IntegerVector group;
        const Rcpp::IntegerVector start;
        const Rcpp::IntegerVector end;
        const Rcpp::NumericMatrix scores;

        std::deque<Range> output;

        inline int size() const { return start.size(); }
        Range get_range(int i) const;

    public:
        IntersectionAggregator(const Rcpp::List arguments);
        void run();
        Rcpp::List get_result();
};

#endif

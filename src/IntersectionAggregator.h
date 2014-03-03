#ifndef PTA_IntersectionAggregator_h
#define PTA_IntersectionAggregator_h
#include <Rcpp.h>
#include <deque>
#include <map>

class IntersectionAggregator {
    private:
        struct Range {
            long group;
            long start;
            long end;
            Rcpp::NumericVector score;

            Range(long g, long s, long e, Rcpp::NumericVector sc) :
                group(g), start(s), end(e), score(sc)
            {}
        };

        const Rcpp::IntegerVector group;
        const Rcpp::IntegerVector start;
        const Rcpp::IntegerVector end;
        const Rcpp::NumericMatrix scores;

        long position;
        std::map<long, Range> open_ranges;
        std::deque<Range> output;

        inline int size() const { return start.size(); }
        inline Rcpp::NumericVector zero_score() { return Rcpp::NumericVector(scores.ncol()); }
        Range get_range(int i) const;
        void close_open_ranges(long until);
        void output_range(long start, long end);

    public:
        IntersectionAggregator(const Rcpp::List arguments);
        void run();
        Rcpp::List get_result();
};

#endif

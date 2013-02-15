#ifndef pta_OptimalPTAProcessor_h
#define pta_OptimalPTAProcessor_h
#include <Rcpp.h>
#include <vector>
#include "vector2.h"
#include "BasePTAProcessor.h"

class OptimalPTAProcessor : public BasePTAProcessor {
    private:
        vector2<double> errors;
        static const int ERRORS_SIZE = 2;

        vector2<int> backindices;

        std::vector<double> cumulative_sums;
        std::vector<double> square_sums;
        std::vector<double> length_sums;

        void merge2(int i);
        void merge_range(int from, int to);
        double sse(int from, int to) const;

    public:
        OptimalPTAProcessor(SEXP start_, SEXP end_, SEXP score_, SEXP count_, SEXP error_, SEXP adjacency_treshold_);
        void run();
};

#endif

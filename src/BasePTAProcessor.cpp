#include "BasePTAProcessor.h"
using namespace Rcpp;

BasePTAProcessor::BasePTAProcessor(SEXP start_, SEXP end_, SEXP score_, SEXP count_, SEXP error_)
    : start(start_), end(end_), score(score_)
{
    count_bound = as<int>(count_);
    error_bound = as<double>(error_);
    if (count_bound) {
        error_bounded = false;
    } else {
        error_bounded = true;

        // calculate maximum error
        double score1 = score[0];
        double length1 = length(0);
        bool first_merge = true;
        int first_i = 0;
        maximum_error = 0;
        minimum_count = 1;
        for (int i = 1; i < size(); ++i) {
            if (adjacent(i - 1, i)) {
                score1 = (length1 * score1 + length(i) * score[i])
                    / (length1 + length(i));
                length1 += length(i);
                maximum_error += pow(length1 * score1 - length(i) * score[i], 2);
                if (first_merge) {
                    maximum_error += pow(length1 * score1 - length(first_i) * score[first_i], 2);
                    first_merge = false;
                }
            } else {
                score1 = score[i];
                length1 = length(i);
                first_merge = true;
                first_i = i;
                ++minimum_count;
            }
        }
    }
}

List BasePTAProcessor::get_result() const {
    return List::create(
            Named("start") = start,
            Named("end") = end,
            Named("score") = score);
}

double BasePTAProcessor::length(int interval) const {
    return end[interval] - start[interval] + 1;
}

bool BasePTAProcessor::adjacent(int i, int j) const {
    return end[i] + 1 == start[j];
}

double BasePTAProcessor::merged_score(int i, int j) const {
    return (length(i) * score[i] + length(j) * score[j])
        / (length(i) + length(j));
}

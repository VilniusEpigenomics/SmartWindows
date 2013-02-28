#include <algorithm>
#include "BasePTAProcessor.h"
using namespace Rcpp;

BasePTAProcessor::BasePTAProcessor(SEXP start_, SEXP end_, SEXP score_, SEXP count_, SEXP error_, SEXP adjacency_treshold_)
    : start(start_), end(end_), score(score_)
{
    count_bound = as<int>(count_);
    error_bound = as<double>(error_);
    adjacency_treshold = as<double>(adjacency_treshold_);
    if (count_bound > 1) {
        error_bounded = false;
    } else {
        error_bounded = true;

        // calculate maximum error
        double score1 = score[0];
        double length1 = length(0);
        int first_i = 0;
        maximum_error = 0;
        minimum_count = 1;
        for (int i = 1; i <= size(); ++i) {
            if (adjacent(i - 1, i) && (i < size())) {
                score1 = (length1 * score1 + length(i) * score[i])
                    / (length1 + length(i));
                length1 += length(i);
            } else {
                for (int j = first_i; j < i; ++j) {
                    maximum_error += length(j) * pow(score1 - score[j], 2);
                }

                if (i < size()) {
                    score1 = score[i];
                    length1 = length(i);
                    first_i = i;
                    ++minimum_count;
                }
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
    if (i == j) {
        return true;
    }

    if (i > j) {
        std::swap(i, j);
    }

    double distance = start[j] - end[i];
    if (distance < 0) {
        stop("Intervals should be sorted and non-overlapping.");
    }
    return distance <= adjacency_treshold;
}

double BasePTAProcessor::merged_score(int i, int j) const {
    return (length(i) * score[i] + length(j) * score[j])
        / (length(i) + length(j));
}

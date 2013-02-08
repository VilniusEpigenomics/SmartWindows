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
        count_bound = size();
        error_bounded = true;
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

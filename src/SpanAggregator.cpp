#include <algorithm>
#include <sstream>
#include <cmath>
#include "SpanAggregator.h"

using namespace Rcpp;

static void assert_is_finite(const NumericVector &x, int row) {
    for (int i = 0; i < x.size(); ++i) {
        double xi = x[i];
        if (!std::isfinite(xi)) {
            std::stringstream msg;
            msg << "Input contains a NaN or infinite value at row " << row + 1 << ", index " << i + 1 << ".";
            throw Rcpp::exception(msg.str().c_str());
        }
    }
}

SpanAggregator::SpanAggregator(const List arguments) :
    start(IntegerVector(static_cast<SEXP>(arguments["start"]))),
    end(IntegerVector(static_cast<SEXP>(arguments["end"]))),
    scores(NumericMatrix(static_cast<SEXP>(arguments["scores"]))),
    span(as<long>(arguments["span"]))
{}

List SpanAggregator::run() {
    std::deque<long> dest_start;
    std::deque<long> dest_end;
    std::deque<NumericVector> dest_scores;

    int i = 0;
    long deststart = start[0];
    long destend = deststart + span;
    NumericVector scores_sum(scores.ncol());
    scores_sum.fill(0);
    long len_sum = 0;
    int count = 0;

    while (i < size()) {
        long intersection_start = std::max(static_cast<long>(start[i]), deststart);
        long intersection_end = std::min(static_cast<long>(end[i]), destend);
        long intersection_len = intersection_end - intersection_start;
        if (intersection_len > 0) {
            // There is an intersection.
            const NumericVector &row = const_cast<NumericMatrix &>(scores)(i, _);
            assert_is_finite(row, i);
            scores_sum += row * static_cast<double>(intersection_len);
            len_sum += intersection_len;
            ++count;
            ++i;
        } else if (end[i] <= deststart) {
            // Current range is left of dest range. Move on.
            ++i;
        } else if (destend <= start[i]) {
            // Current range if right of dest range. Append resulting span if any
            // ranges were merged.
            if (count > 0) {
                dest_start.push_back(deststart);
                dest_end.push_back(destend);
                dest_scores.push_back(scores_sum / static_cast<double>(len_sum));
                scores_sum.fill(0);
                len_sum = 0;
                count = 0;
            }
            deststart += span;
            destend += span;
        }
    }

    if (count > 0) {
        dest_start.push_back(deststart);
        dest_end.push_back(destend);
        dest_scores.push_back(scores_sum / static_cast<double>(len_sum));
    }

    long result_size = dest_start.size();

    NumericVector result_start(result_size);
    NumericVector result_end(result_size);
    NumericMatrix result_scores(result_size, scores.ncol());

    std::deque<long>::const_iterator start_iter = dest_start.begin();
    std::deque<long>::const_iterator end_iter = dest_end.begin();
    std::deque<NumericVector>::const_iterator scores_iter = dest_scores.begin();
    for (int result_i = 0; result_i < result_size; ++result_i) {
        result_start[result_i] = *start_iter;
        result_end[result_i] = *end_iter;
        result_scores(result_i, _) = *scores_iter;
        ++start_iter;
        ++end_iter;
        ++scores_iter;
    }

    List result = List::create(
            Named("start") = result_start,
            Named("end") = result_end,
            Named("scores") = result_scores);
    return result;
}

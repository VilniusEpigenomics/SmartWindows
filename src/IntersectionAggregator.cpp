#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include "IntersectionAggregator.h"

using namespace Rcpp;

RcppExport SEXP intersectionAggregate(SEXP arguments_) {
BEGIN_RCPP
    const List arguments(arguments_);
    IntersectionAggregator a(arguments);
    a.run();
    return a.get_result();
END_RCPP
}

IntersectionAggregator::IntersectionAggregator(const List arguments) :
    group(IntegerVector(static_cast<SEXP>(arguments["group"]))),
    start(IntegerVector(static_cast<SEXP>(arguments["start"]))),
    end(IntegerVector(static_cast<SEXP>(arguments["end"]))),
    scores(NumericMatrix(static_cast<SEXP>(arguments["scores"]))),
    position(std::numeric_limits<long>::min())
{}

IntersectionAggregator::Range IntersectionAggregator::get_range(int i) const {
    const NumericVector s = const_cast<NumericMatrix&>(scores)(i, _);
    return Range(group[i], start[i], end[i], s);
}

void IntersectionAggregator::close_open_ranges(long until) {
    while (!open_ranges.empty() && (position < until)) {
        long first_open_end = open_ranges.begin()->first;
        long r_end = std::min(until, first_open_end);
        output_range(position, r_end);
        if (first_open_end <= until) {
            open_ranges.erase(first_open_end);
        }
        position = r_end;
    }
}

void IntersectionAggregator::output_range(long start, long end) {
    NumericVector sum = zero_score();
    int count = 0;
    long group;
    for (std::map<long, Range>::const_iterator it = open_ranges.begin(); it != open_ranges.end(); ++it) {
        const Range &r = it->second;
        group = r.group;
        sum += r.score;
        count++;
    }
    const NumericVector score = sum / static_cast<double>(count);
    output.push_back(Range(group, start, end, score));
}

void IntersectionAggregator::run() {
    for (int i = 0; i < size(); ++i) {
        const Range r = get_range(i);
        close_open_ranges(r.start);
        open_ranges.insert(std::make_pair(r.end, r));
        position = r.start;
    }
    close_open_ranges(std::numeric_limits<long>::max());
}

List IntersectionAggregator::get_result() {
    int output_size = output.size();
    IntegerVector output_group(output_size);
    IntegerVector output_start(output_size);
    IntegerVector output_end(output_size);
    NumericMatrix output_scores(output_size, scores.ncol());
    std::deque<Range>::const_iterator output_iter = output.begin();
    for (int output_i = 0; output_i < output_size; ++output_i) {
        const Range &r = *output_iter;
        output_group[output_i] = r.group;
        output_start[output_i] = r.start;
        output_end[output_i] = r.end;
        output_scores(output_i, _) = r.score;
        ++output_iter;
    }

    List result = List::create(
            Named("group") = output_group,
            Named("start") = output_start,
            Named("end") = output_end,
            Named("scores") = output_scores);
    return result;
}

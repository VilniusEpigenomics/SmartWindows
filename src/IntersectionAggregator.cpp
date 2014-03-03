#include <algorithm>
#include <sstream>
#include <cassert>
#include <cmath>
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
    scores(NumericMatrix(static_cast<SEXP>(arguments["scores"])))
{}

IntersectionAggregator::Range IntersectionAggregator::get_range(int i) const {
    Range r;
    r.group = group[i];
    r.start = start[i];
    r.end = end[i];
    r.score = scores[i];
    return r;
}

void IntersectionAggregator::run() {
    for (int i = 0; i < size(); ++i) {
        const Range r = get_range(i);
        output.push_back(r);
    }
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

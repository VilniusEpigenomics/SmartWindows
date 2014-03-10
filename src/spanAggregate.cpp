#include <deque>
#include <sstream>
#include <Rcpp.h>

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

RcppExport SEXP spanAggregate(SEXP arguments_) {
BEGIN_RCPP
    const List arguments(arguments_);
    const Rcpp::IntegerVector start(IntegerVector(static_cast<SEXP>(arguments["start"])));
    const Rcpp::IntegerVector end(IntegerVector(static_cast<SEXP>(arguments["end"])));
    const Rcpp::NumericMatrix scores(NumericMatrix(static_cast<SEXP>(arguments["scores"])));
    const long span = as<long>(arguments["span"]);
    const int ranges_count = start.size();

    std::deque<long> dest_start;
    std::deque<long> dest_end;
    std::deque< std::vector<double> > dest_scores;

    int i = 0;
    long spanstart = start[0];
    long spanend = spanstart + span;
    NumericVector scores_sum(scores.ncol());
    scores_sum.fill(0);
    long len_sum = 0;
    int count = 0;

    long whole_start = start[0];
    long whole_end = end[end.size() - 1];

    while (spanend < whole_end + span) {
        while (spanend <= start[i]) {
            // Current range if right of span. Move to next span.
            spanstart += span;
            spanend += span;
        }

        while (end[i] <= spanstart) {
            // Current range is left of span. Move to next range.
            ++i;
        }

merge_range:
        {
            // There is an intersection.
            long intersection_start = std::max(static_cast<long>(start[i]), spanstart);
            long intersection_end = std::min(static_cast<long>(end[i]), spanend);
            long intersection_len = intersection_end - intersection_start;
            const NumericVector &row = const_cast<NumericMatrix &>(scores)(i, _);
            assert_is_finite(row, i);
            scores_sum += row * static_cast<double>(intersection_len);
            len_sum += intersection_len;
            ++count;

            if ((end[i] < spanend) && (i < ranges_count - 1)) {
                ++i;
                goto merge_range;
            }
        }

        if (count > 0) {
            dest_start.push_back(spanstart);
            dest_end.push_back(spanend);
            const NumericVector mean_score = scores_sum / static_cast<double>(len_sum);
            dest_scores.push_back(as< std::vector<double> >(mean_score));
            scores_sum.fill(0);
            len_sum = 0;
            count = 0;
        }
        spanstart += span;
        spanend += span;
    }

    long result_size = dest_start.size();

    NumericVector result_start(result_size);
    NumericVector result_end(result_size);
    NumericMatrix result_scores(result_size, scores.ncol());

    std::deque<long>::const_iterator start_iter = dest_start.begin();
    std::deque<long>::const_iterator end_iter = dest_end.begin();
    std::deque< std::vector<double> >::const_iterator scores_iter = dest_scores.begin();
    for (int result_i = 0; result_i < result_size; ++result_i) {
        result_start[result_i] = *start_iter;
        result_end[result_i] = *end_iter;
        const NumericVector score(scores_iter->begin(), scores_iter->end());
        result_scores(result_i, _) = score;
        ++start_iter;
        ++end_iter;
        ++scores_iter;
    }

    List result = List::create(
            Named("start") = result_start,
            Named("end") = result_end,
            Named("scores") = result_scores);
    return result;
END_RCPP
}

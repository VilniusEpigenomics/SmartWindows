#include "OptimalPTAProcessor.h"
using namespace Rcpp;
using namespace std;

OptimalPTAProcessor::OptimalPTAProcessor(SEXP start_, SEXP end_, SEXP score_, SEXP count_, SEXP error_) :
    BasePTAProcessor(start_, end_, score_, count_, error_)
{
    minimum_count = 1;
    for (int i = 0; i < size() - 1; ++i) {
        if (!adjacent(i, i + 1)) {
            ++minimum_count;
            nonadjacencies.push_back(i);
        }
    }

    cumulative_sums.resize(size());
    cumulative_sums[0] = length(0) * score[0];
    for (int i = 1; i < size(); ++i) {
        cumulative_sums[i] = cumulative_sums[i - 1] + length(i) * score[i];
    }

    square_sums.resize(size());
    square_sums[0] = length(0) * pow(score[0], 2);
    for (int i = 1; i < size(); ++i) {
        square_sums[i] = square_sums[i - 1] + length(i) * pow(score[i], 2);
    }

    length_sums.resize(size());
    length_sums[0] = length(0);
    for (int i = 1; i < size(); ++i) {
        length_sums[i] = length_sums[i - 1] + length(i);
    }

    errors.assign(2, size(), INFINITY);
    maximum_error = sse(minimum_count, size() - 1);

    backindices.assign(count_bound, size(), -1);
}

void OptimalPTAProcessor::merge2(int i) {
    score[i] = merged_score(i, i + 1);
    end[i] = end[i + 1];

    start.erase(i + 1);
    end.erase(i + 1);
    score.erase(i + 1);
}

void OptimalPTAProcessor::merge_range(int from, int to) {
    while (to > from) {
        merge2(to - 1);
        --to;
    }
}

double OptimalPTAProcessor::sse(int from, int to) const {
    double s_sub, ss_sub, l_sub;
    if (from > 0) {
        s_sub = cumulative_sums[from - 1];
        ss_sub = square_sums[from - 1];
        l_sub = length_sums[from - 1];
    } else {
        s_sub = 0;
        ss_sub = 0;
        l_sub = 0;
    }
    return square_sums[to] - ss_sub - pow(cumulative_sums[to] - s_sub, 2)/(length_sums[to] - l_sub);
}

void OptimalPTAProcessor::run() {
    int c = count_bound;
    for (int k = 0; k < c; ++k) {
        for (int i = k; i < size(); ++i) {
            if (k == 0) {
                errors(0, i) = sse(0, i);
            } else {
                for (int j = i - 1; j >= k - 1; --j) {
                    double err1 = errors((k - 1) % ERRORS_SIZE, j);
                    double err2 = sse(j, i);
                    if (err1 + err2 < errors(k % ERRORS_SIZE, i)) {
                        errors(k % ERRORS_SIZE, i) = err1 + err2;
                        backindices(k, i) = j;
                    }
                    if (err2 > errors(k % ERRORS_SIZE, i)) break;
                }
            }
        }

        if (error_bounded && (errors(k % ERRORS_SIZE, size() - 1) <= error_bound * maximum_error)) {
            c = k;
            break;
        }
    }

    int n = size() - 1;
    while ((c > 0) && (n > 0)) {
        int j = backindices(c - 1, n);
        merge_range(j + 1, n);
        n = j;
        --c;
    }
}

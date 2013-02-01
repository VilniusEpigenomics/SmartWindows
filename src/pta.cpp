#include "pta.h"
#include <iostream>
#include <cmath>
#include <limits>
#include <vector>
#include <queue>
using namespace Rcpp;
using namespace std;

template<typename T>
struct vector2 {
    int rows;
    int columns;
    vector<T> data;

    void assign(int rows_, int columns_, T value) {
        rows = rows_;
        columns = columns_;
        data.assign(rows * columns, value);
    }

    T operator()(int i, int j) const {
        return data[i * columns + j];
    }

    T& operator()(int i, int j) {
        return data[i * columns + j];
    }
};

struct BasePTAProcessor {
    NumericVector start;
    NumericVector end;
    NumericVector score;

    int count_bound;
    double error_bound;
    bool error_bounded;

    BasePTAProcessor(SEXP start_, SEXP end_, SEXP score_, SEXP count_, SEXP error_) :
        start(start_), end(end_), score(score_)
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

    List get_result() const {
        return List::create(
                Named("start") = start,
                Named("end") = end,
                Named("score") = score);
    }

    inline int size() const {
        return start.size();
    }

    double length(int interval) const {
        return end[interval] - start[interval] + 1;
    }

    bool adjacent(int i, int j) const {
        return end[i] + 1 == start[j];
    }

    double merged_score(int i, int j) const {
        return (length(i) * score[i] + length(j) * score[j])
             / (length(i) + length(j));
    }
};

struct OptimalPTAProcessor : BasePTAProcessor {
    int minimum_count;
    vector<int> nonadjacencies;

    vector2<double> errors;
    static const int ERRORS_SIZE = 2;
    double maximum_error;

    vector2<int> backindices;

    vector<double> cumulative_sums;
    vector<double> square_sums;
    vector<double> length_sums;

    OptimalPTAProcessor(SEXP start_, SEXP end_, SEXP score_, SEXP count_, SEXP error_) :
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

    void merge2(int i) {
        score[i] = merged_score(i, i + 1);
        end[i] = end[i + 1];

        start.erase(i + 1);
        end.erase(i + 1);
        score.erase(i + 1);
    }

    void merge_range(int from, int to) {
        while (to > from) {
            merge2(to - 1);
            --to;
        }
    }

    double sse(int from, int to) const {
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

    void run() {
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
};

struct GreedyPTAProcessor : BasePTAProcessor {
    struct SimilarityComparator {
        const GreedyPTAProcessor& processor;

        SimilarityComparator(const GreedyPTAProcessor& processor_) :
            processor(processor_) {}

        double operator()(int i, int j) {
            return processor.key(i) > processor.key(j);
        }
    };

    SimilarityComparator comparator;
    priority_queue<int, vector<int>, SimilarityComparator> heap;
    int node_count;
    int first_node;
    int last_node;
    vector<int> next_node;
    vector<int> prev_node;

    void insert(int node) {
        heap.push(node);
    }

    int peek() const {
        return heap.top();
    }

    double key(int node) const {
        if (node == first_node) {
            return INFINITY;
        } else {
            return dsim(prev_node[node], node);
        }
    }

    void erase(int node) {
        int prev = prev_node[node];
        int next = next_node[node];

        if (node == first_node) {
            first_node = next;
        } else {
            next_node[prev] = next;
        }

        if (node == last_node) {
            last_node = prev;
        } else {
            prev_node[next] = prev;
        }

        --node_count;
    }

    bool merge() {
        if (heap.size() <= 1) return false;
        int top = peek();
        if (key(top) == INFINITY) return false;

        int prev = prev_node[top];
        score[prev] = merged_score(prev, top);
        erase(top);
        heap.pop();
        resort();
        return true;
    }

    void resort() {
        make_heap(
                const_cast<int*>(&heap.top()),
                const_cast<int*>(&heap.top()) + heap.size(),
                comparator);
    }

    GreedyPTAProcessor(SEXP start_, SEXP end_, SEXP score_, SEXP count_, SEXP error_) :
        BasePTAProcessor(start_, end_, score_, count_, error_),
        comparator(SimilarityComparator(*this)),
        heap(comparator)
    {
        node_count = size();
        prev_node.resize(node_count);
        next_node.resize(node_count);
        first_node = 0;
        last_node = node_count - 1;
        for (int i = 0; i < size(); ++i) {
            prev_node[i] = i - 1;
            next_node[i] = i + 1;
            insert(i);
        }
        next_node[last_node] = -1;
    }

    double dsim(int i, int j) const {
        if (adjacent(i, j)) {
            const double z = merged_score(i, j);
            return length(i) * pow(z - score[i], 2) + length(j) * pow(z - score[j], 2);
        } else {
            return INFINITY;
        }
    }

    void run() {
        cout << "heap size: " << heap.size() << endl;

        for (int i = 0; i < size() - 1; ++i) {
            if (!adjacent(i, i + 1)) {
                cout << "nonadjacency at " << i << endl;
            }
        }

        while (node_count > count_bound) {
            if (!merge()) {
                cout << "merge broken with " << node_count << " nodes" << endl;
                break;
            }
        }

        NumericVector newstart(node_count);
        NumericVector newend(node_count);
        NumericVector newscore(node_count);
        int node = first_node;
        int i = 0;
        while (node >= 0) {
            newstart[i] = start[node];
            newend[i] = end[node];
            newscore[i] = score[node];
            node = next_node[node];
            ++i;
        }

        start = newstart;
        end = newend;
        score = newscore;
    }
};

SEXP PTA(SEXP start, SEXP end, SEXP score, SEXP count, SEXP error) {
    OptimalPTAProcessor p(start, end, score, count, error);
    p.run();
    return p.get_result();
}

SEXP gPTA(SEXP start, SEXP end, SEXP score, SEXP count, SEXP error) {
    GreedyPTAProcessor p(start, end, score, count, error);
    p.run();
    return p.get_result();
}

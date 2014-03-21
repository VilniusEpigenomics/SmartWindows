#include <algorithm>
#include <sstream>
#include <cassert>
#include <cmath>
#include "PTAProcessor.h"

using namespace Rcpp;

#define FOR_EACH_HEAP(var) for (int (var) = 0; (var) < nheaps; ++(var))

#ifdef NDEBUG
#define CHECK_HEAP
#else
#define CHECK_HEAP \
    for (int __CHECK_HEAP_heap = 0; __CHECK_HEAP_heap < nheaps; ++__CHECK_HEAP_heap) \
        assert(std::is_heap(heaps[__CHECK_HEAP_heap].begin(), heaps[__CHECK_HEAP_heap].end(), greaters[__CHECK_HEAP_heap]));
#endif

/// UTILITY FUNCTIONS

static inline double sum_sq(const NumericVector& x) {
    double sum = 0.0;
    for (int i = 0; i < x.size(); ++i) {
        double x_i = x[i];
        sum += x_i * x_i;
    }
    return sum;
}

struct LessThanIndirect {
    const NumericVector& x;
    LessThanIndirect(const NumericVector& x_) : x(x_) {}
    double operator()(int a, int b) const {
        return x[a] < x[b];
    }
};

inline static NumericVector rank(const NumericVector& x) {
    NumericVector order(x.size());
    for (int i = 0; i < order.size(); ++i) {
        order[i] = i;
    }

    const LessThanIndirect less_than(x);
    std::sort(order.begin(), order.end(), less_than);

    NumericVector r(x.size());
    for (int i = 0; i < r.size(); ++i) {
        r[order[i]] = i;
    }
    return r;
}

static inline double correlation(const NumericVector &x, const NumericVector &y, bool spearman = false) {
    const int n = x.size();

    if (spearman) {
        // Ties in ranks are improbable (and not handled by the ranking
        // algorithm), so use simplified method to calculate rho.
        const NumericVector diff = rank(x) - rank(y);
        return 1 - (6 * sum_sq(diff)) / (n * (n * n - 1));
    } else {
        const double mean_x = mean(x);
        const double mean_y = mean(y);

        double var_x = 0;
        double var_y = 0;
        double covar = 0;

        for (int i = 0; i < n; ++i) {
            const double delta_x = x[i] - mean_x;
            const double delta_y = y[i] - mean_y;
            var_x += delta_x * delta_x;
            var_y += delta_y * delta_y;
            covar += delta_x * delta_y;
        }

        const double max_covar = std::sqrt(var_x * var_y);

        if (max_covar == 0) {
            return 0;
        }

        return covar / max_covar;
    }
}

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


/// METHODS

double PTAProcessor::length(int interval) const {
    return end[interval] - start[interval] + 1;
}

bool PTAProcessor::adjacent(int i, int j) const {
    if (i == j) {
        return true;
    }

    if (i > j) {
        std::swap(i, j);
    }

    const double distance = start[j] - end[i];
    if (distance < 0) {
        throw Rcpp::exception("Intervals should be sorted and non-overlapping.");
    }
    return distance <= adjacency_threshold;
}

NumericVector PTAProcessor::merged_scores(int i, int j) const {
    const NumericVector scores_i = const_cast<PTAProcessor*>(this)->scores(i, _);
    const NumericVector scores_j = const_cast<PTAProcessor*>(this)->scores(j, _);
    return (length(i) * scores_i + length(j) * scores_j)
        / (length(i) + length(j));
}

double PTAProcessor::key(int heap, int nodeid) const {
    int previd = nodeid;
    for (int i = -1; i < heap; ++i) {
        if (previd == first_node) return INFINITY;
        previd = nodes[previd].prev;
    }

    if (!adjacent(previd, nodeid)) return INFINITY;

    if (correlation_mode) {
        double cor = node_correlation(previd, nodeid);
        return 1 - cor;
    } else {
        return merge_error(previd, nodeid);
    }
}

int PTAProcessor::parent(int i) const {
    return (i - 1)/2;
}

int PTAProcessor::left_child(int i) const {
    return 2 * i + 1;
}

int PTAProcessor::right_child(int i) const {
    return 2 * i + 2;
}

void PTAProcessor::heap_swap(int heap, int i, int j) {
    std::swap(heaps[heap][i], heaps[heap][j]);
    nodes[heaps[heap][i]].positions[heap] = i;
    nodes[heaps[heap][j]].positions[heap] = j;
}

void PTAProcessor::heap_up(int heap, int i) {
    while (i > 0) {
        int par = parent(i);
        if (greaters[heap](heaps[heap][par], heaps[heap][i])) {
            heap_swap(heap, i, par);
            i = par;
        } else {
            break;
        }
    }
}

void PTAProcessor::heap_down(int heap, int i) {
    while (true) {
        int smallest = i;
        int right = right_child(i);
        if (right < static_cast<int>(heaps[heap].size()) && greaters[heap](heaps[heap][smallest], heaps[heap][right])) {
            smallest = right;
        }
        int left = left_child(i);
        if (left < static_cast<int>(heaps[heap].size()) && greaters[heap](heaps[heap][smallest], heaps[heap][left])) {
            smallest = left;
        }

        if (smallest == i) {
            break;
        } else {
            heap_swap(heap, i, smallest);
            i = smallest;
        }
    }
}

void PTAProcessor::heap_delete(int heap, int i) {
    heaps[heap][i] = heaps[heap].back();
    nodes[heaps[heap][i]].positions[heap] = i;
    heaps[heap].pop_back();

    if (i > 0 && greaters[heap](heaps[heap][parent(i)], heaps[heap][i])) {
        heap_up(heap, i);
    } else {
        heap_down(heap, i);
    }
}

void PTAProcessor::heap_insert(int heap, int nodeid) {
    heaps[heap].push_back(nodeid);
    int i = heaps[heap].size() - 1;
    nodes[nodeid].positions[heap] = i;
    heap_up(heap, i);
}

void PTAProcessor::update_node(int nodeid) {
    if (nodeid == -1) return;
    Node& node = nodes[nodeid];
    FOR_EACH_HEAP(heap) {
        heap_delete(heap, node.positions[heap]);
        node.keys[heap] = key(heap, node.id);
        heap_insert(heap, node.id);
    }
}

bool PTAProcessor::merge(int minheap, int minnode, double *error) {
    const Node& top = nodes[minnode];
    if (top.keys[minheap] == INFINITY) return false;

    int nodeid = minnode;
    for (int i = -1; i < minheap; ++i) {
        Node& node = nodes[nodeid];
        FOR_EACH_HEAP(j) {
            heap_delete(j, node.positions[j]);
        }
        node.alive = false;
        if (i >= 0) {
            node.skipped = true;
        }
        if (nodeid == first_node) return false;
        nodeid = nodes[nodeid].prev;
    }
    CHECK_HEAP;
    node_count -= minheap + 1;
    Node& repl = nodes[nodeid]; // replacement

    const NumericVector merged = merged_scores(repl.id, minnode);
    if (error) {
        int i = repl.id;
        int j = minnode;
        const NumericVector diff_i = merged - const_cast<PTAProcessor *>(this)->scores(i, _);
        const NumericVector diff_j = merged - const_cast<PTAProcessor *>(this)->scores(j, _);
        *error = length(i) * sum_sq(diff_i) + length(j) * sum_sq(diff_j);
    }

    scores(repl.id, _) = merged;
    end[repl.id] = end[minnode];
    repl.next = top.next;
    if (repl.next != -1) nodes[repl.next].prev = repl.id;

    nodeid = repl.id;
    for (int i = -1; i < nheaps; ++i) {
        Node& node = nodes[nodeid];
        update_node(node.id);
        if (node.next == -1) break;
        nodeid = node.next;
    }
    CHECK_HEAP;

    return true;
}

PTAProcessor::PTAProcessor(const List arguments) :
    original_start(NumericVector(static_cast<SEXP>(arguments["start"]))),
    original_end(NumericVector(static_cast<SEXP>(arguments["end"]))),
    original_scores(NumericMatrix(static_cast<SEXP>(arguments["scores"]))),
    count_bound(as<int>(arguments["countBound"])),
    error_bound(as<double>(arguments["errorBound"])),
    cumulative_error_bound(as<double>(arguments["cumulativeErrorBound"])),
    adjacency_threshold(as<double>(arguments["adjacencyThreshold"])),
    mode(as<int>(arguments["mode"])),
    correlation_mode(as<int>(arguments["mode"]) != PTA_MODE_NORMAL),
    correlation_spearman(as<int>(arguments["mode"]) == PTA_MODE_CORRELATION_SPEARMAN),
    correlation_bound(as<double>(arguments["correlationBound"])),
    nheaps(as<int>(arguments["skip"]) + 1)
{
    start = clone(original_start);
    end = clone(original_end);
    scores = clone(original_scores);

    switch (mode) {
        case PTA_MODE_NORMAL:
        case PTA_MODE_CORRELATION:
        case PTA_MODE_CORRELATION_SPEARMAN:
            break;
        default:
            throw Rcpp::exception("Bad PTA mode.");
    }

    if (cumulative_error_bound == INFINITY) {
        maximum_error = INFINITY;
    } else {
        // calculate maximum error
        NumericVector score1 = clone(NumericVector(scores(0, _)));
        double length1 = length(0);
        int first_i = 0;
        maximum_error = 0;
        minimum_count = 1;
        for (int i = 1; i <= size(); ++i) {
            if ((i < size()) && adjacent(i - 1, i)) {
                score1 = (length1 * score1 + length(i) * scores(i, _))
                       / (length1 + length(i));
                length1 += length(i);
            } else {
                for (int j = first_i; j < i; ++j) {
                    const NumericVector diff = score1 - scores(j, _);
                    maximum_error += length(j) * sum_sq(diff);
                }

                if (i < size()) {
                    score1 = scores(i, _);
                    length1 = length(i);
                    first_i = i;
                    ++minimum_count;
                }
            }
        }
    }

    node_count = size();
    nodes.resize(node_count);
    heaps.resize(nheaps);
    FOR_EACH_HEAP(i) heaps[i].resize(node_count);
    first_node = 0;
    last_node = node_count - 1;
    for (int i = 0; i < node_count; ++i) {
        assert_is_finite(scores(i, _), i);
        Node& node = nodes[i];
        node.alive = true;
        node.skipped = false;
        node.prev = i - 1;
        node.id = i;
        node.next = i + 1;
        node.keys.resize(nheaps);
        node.positions.resize(nheaps);
        FOR_EACH_HEAP(heap) {
            node.keys[heap] = key(heap, node.id);
            heaps[heap][i] = i;
        }
    }
    nodes[last_node].next = -1;
    greaters.resize(nheaps);
    FOR_EACH_HEAP(heap) {
        greaters[heap].set(this, heap);
        std::make_heap(heaps[heap].begin(), heaps[heap].end(), greaters[heap]);
        for (int i = 0; i < node_count; ++i) {
            nodes[heaps[heap][i]].positions[heap] = i;
        }
    }
}

double PTAProcessor::merge_error(int i, int j) const {
    const NumericVector z = merged_scores(i, j);
    const NumericVector diff_i = z - const_cast<PTAProcessor *>(this)->scores(i, _);
    const NumericVector diff_j = z - const_cast<PTAProcessor *>(this)->scores(j, _);
    return length(i) * sum_sq(diff_i) + length(j) * sum_sq(diff_j);
}

double PTAProcessor::node_correlation(int x, int y) const {
    const NumericVector scores_x = const_cast<PTAProcessor*>(this)->scores(x, _);
    const NumericVector scores_y = const_cast<PTAProcessor*>(this)->scores(y, _);
    return correlation(scores_x, scores_y, correlation_spearman);
}

List PTAProcessor::run() {
    const double abs_error_bound = cumulative_error_bound * maximum_error;
    double cumulative_error = 0;
    while (node_count > count_bound) {
        int minheap = 0;
        int minnode = heaps[0].front();
        double minkey = key(minheap, minnode);
        for (int i = 1; i < nheaps; ++i) {
            if (key(i, heaps[i].front()) < minkey) {
                minheap = i;
                minnode = heaps[i].front();
                minkey = key(minheap, minnode);
            }
        }

        if (correlation_mode) {
           if (1 - minkey <= correlation_bound) {
               break;
           }
        }

        double error;
        if (!merge(minheap, minnode, &error)) {
            break;
        }

        cumulative_error += error;

        if ((error > error_bound) ||
            (cumulative_error > abs_error_bound))
        {
            break;
        }
    }

    NumericVector newstart(node_count);
    NumericVector newend(node_count);
    NumericMatrix newscores(node_count, scores.ncol());
    IntegerVector groups(original_start.size());
    int groupid = -1;
    for (int i = 0; i < original_start.size(); ++i) {
        if (nodes[i].alive) {
            ++groupid;
            newstart[groupid] = start[i];
            newend[groupid] = end[i];
            newscores(groupid, _) = scores(i, _);
        }
        if (nodes[i].skipped) {
            groups[i] = -1;
        } else {
            groups[i] = groupid;
        }
    }

    List result = List::create(
            Named("start") = newstart,
            Named("end") = newend,
            Named("scores") = newscores,
            Named("groups") = groups,
            Named("cumulativeError") = cumulative_error);
    return result;
}

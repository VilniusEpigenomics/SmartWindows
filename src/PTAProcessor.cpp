#include <algorithm>
#include <cassert>
#include "PTAProcessor.h"
using namespace Rcpp;
using namespace std;

#define FOR_EACH_HEAP(var) for (int (var) = 0; (var) < nheaps; ++(var))

#ifdef NDEBUG
#define CHECK_HEAP
#else
#define CHECK_HEAP \
    for (int __CHECK_HEAP_heap = 0; __CHECK_HEAP_heap < nheaps; ++__CHECK_HEAP_heap) \
        assert(is_heap(heaps[__CHECK_HEAP_heap].begin(), heaps[__CHECK_HEAP_heap].end(), greaters[__CHECK_HEAP_heap]));
#endif

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

    double distance = start[j] - end[i];
    if (distance < 0) {
        throw Rcpp::exception("Intervals should be sorted and non-overlapping.");
    }
    return distance <= adjacency_threshold;
}

template <typename T>
struct AffineTransform {
    double intercept, coefficient;
    AffineTransform(double a, double b) : intercept(a), coefficient(b) {}

    inline T operator()(const T& x) const {
        return intercept + coefficient * x;
    }

    inline AffineTransform<T> inverse() const {
        return AffineTransform<T>(-(intercept / coefficient), 1 / coefficient);
    }
};

static inline AffineTransform<NumericVector> linear_regression(const NumericVector& x, const NumericVector& y) {
    double mean_x = mean(x);
    double mean_y = mean(y);
    double b = (mean(x * y) - mean_x * mean_y) / (mean(x*x) - mean_x * mean_x);
    double a = mean_y - b * mean_x;
    return AffineTransform<NumericVector>(a, b);
}

NumericVector PTAProcessor::merged_scores(int i, int j) const {
    const NumericVector scores_i = const_cast<PTAProcessor*>(this)->scores(i, _);
    const NumericVector scores_j = const_cast<PTAProcessor*>(this)->scores(j, _);
    if (!((mode == PTA_MODE_CORRELATION) && correlation_newmerge)) {
        return (length(i) * scores_i + length(j) * scores_j)
             / (length(i) + length(j));
    } else {
        NumericVector scores_i_m = clone(scores_i);
        NumericVector scores_j_m = clone(scores_j);
        if (length(i) >= length(j)) {
            scores_j_m = linear_regression(scores_i_m, scores_j_m).inverse()(scores_j_m);
        } else {
            scores_i_m = linear_regression(scores_j_m, scores_i_m).inverse()(scores_i_m);
        }
        return (length(i) * scores_i_m + length(j) * scores_j_m)
             / (length(i) + length(j));
    }
}

double PTAProcessor::key(int heap, int nodeid) const {
    int previd = nodeid;
    for (int i = -1; i < heap; ++i) {
        if (previd == first_node) return INFINITY;
        previd = nodes[previd].prev;
    }

    if (!adjacent(previd, nodeid)) return INFINITY;

    switch (mode) {
        case PTA_MODE_NORMAL:
            return dsim(previd, nodeid);
        case PTA_MODE_CORRELATION:
            {
                double cor = correlation(previd, nodeid);
                if (correlation_newmerge) {
                    return 1 - abs(cor);
                } else {
                    return 1 - cor;
                }
            }
        default:
            throw Rcpp::exception("Bad PTA mode.");
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
    swap(heaps[heap][i], heaps[heap][j]);
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

bool PTAProcessor::merge(int minheap, int minnode) {
    const Node& top = nodes[minnode];
    if (top.keys[minheap] == INFINITY) return false;

    int nodeid = minnode;
    for (int i = -1; i < minheap; ++i) {
        Node& node = nodes[nodeid];
        FOR_EACH_HEAP(j) {
            heap_delete(j, node.positions[j]);
        }
        if (nodeid == first_node) return false;
        nodeid = nodes[nodeid].prev;
    }
    CHECK_HEAP;
    node_count -= minheap + 1;
    Node& repl = nodes[nodeid]; // replacement

    scores(repl.id, _) = merged_scores(repl.id, minnode);
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

static inline double sum_sq(const NumericVector& x) {
    double sum = 0.0;
    for (int i = 0; i < x.size(); ++i) {
        double x_i = x[i];
        sum += x_i * x_i;
    }
    return sum;
}

PTAProcessor::PTAProcessor(
        SEXP start_,
        SEXP end_,
        SEXP scores_,
        SEXP count_bound_,
        SEXP error_bound_,
        SEXP adjacency_threshold_,
        SEXP skip_,
        SEXP mode_,
        SEXP correlation_bound_,
        SEXP correlation_spearman_,
        SEXP correlation_newmerge_)
    : start(clone(start_)), end(clone(end_)), scores(clone(scores_))
{
    mode = as<int>(mode_);

    count_bound = as<int>(count_bound_);
    error_bound = as<double>(error_bound_);
    correlation_bound = as<double>(correlation_bound_);
    correlation_spearman = as<int>(correlation_spearman_);
    correlation_newmerge = as<int>(correlation_newmerge_);
    adjacency_threshold = as<double>(adjacency_threshold_);
    if ((count_bound > 1) || (mode == PTA_MODE_CORRELATION)) {
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

    nheaps = as<int>(skip_) + 1;
    node_count = size();
    nodes.resize(node_count);
    heaps.resize(nheaps);
    FOR_EACH_HEAP(i) heaps[i].resize(node_count);
    first_node = 0;
    last_node = node_count - 1;
    for (int i = 0; i < node_count; ++i) {
        Node& node = nodes[i];
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
        make_heap(heaps[heap].begin(), heaps[heap].end(), greaters[heap]);
        for (int i = 0; i < node_count; ++i) {
            nodes[heaps[heap][i]].positions[heap] = i;
        }
    }
}

double PTAProcessor::dsim(int i, int j) const {
    const NumericVector z = merged_scores(i, j);
    const NumericVector diff_i = z - const_cast<PTAProcessor *>(this)->scores(i, _);
    const NumericVector diff_j = z - const_cast<PTAProcessor *>(this)->scores(j, _);
    return length(i) * sum_sq(diff_i) + length(j) * sum_sq(diff_j);
}

struct LessThanIndirect {
    const NumericVector& x;
    LessThanIndirect(const NumericVector& x_) : x(x_) {}
    double operator()(int a, int b) const {
        return x[a] < x[b];
    }
};

static NumericVector rank(const NumericVector& x) {
    NumericVector r(x.size());
    for (int i = 0; i < r.size(); ++i) {
        r[i] = i;
    }

    const LessThanIndirect less_than(x);
    sort(r.begin(), r.end(), less_than);
    return r;
}

double PTAProcessor::correlation(int x, int y) const {
    NumericVector scores_x;
    NumericVector scores_y;

    if (correlation_spearman) {
        scores_x = rank(const_cast<PTAProcessor *>(this)->scores(x, _));
        scores_y = rank(const_cast<PTAProcessor *>(this)->scores(y, _));
    } else {
        scores_x = const_cast<PTAProcessor *>(this)->scores(x, _);
        scores_x = const_cast<PTAProcessor *>(this)->scores(y, _);
    }

    int n = scores_x.size();

    double mean_x = mean(scores_x);
    double mean_y = mean(scores_y);

    double var_x = 0;
    double var_y = 0;
    double covar = 0;

    for (int i = 0; i < n; ++i) {
        double delta_x = scores_x[i] - mean_x;
        double delta_y = scores_y[i] - mean_y;
        var_x += delta_x * delta_x;
        var_y += delta_y * delta_y;
        covar += delta_x * delta_y;
    }

    return covar / sqrt(var_x * var_y);
}

List PTAProcessor::run() {
    double abs_error_bound = error_bound * maximum_error;
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

        if (mode == PTA_MODE_CORRELATION) {
           if (1 - minkey <= correlation_bound) {
               break;
           }
        } else {
            if (cumulative_error + minkey > abs_error_bound) {
                break;
            }
        }

        if (!merge(minheap, minnode)) {
            break;
        }

        cumulative_error += minkey;
    }

    NumericVector newstart(node_count);
    NumericVector newend(node_count);
    NumericMatrix newscores(node_count, scores.ncol());
    int nodeid = first_node;
    int i = 0;
    while (nodeid >= 0) {
        newstart[i] = start[nodeid];
        newend[i] = end[nodeid];
        newscores(i, _) = scores(nodeid, _);
        nodeid = nodes[nodeid].next;
        ++i;
    }

    return List::create(
            Named("start") = newstart,
            Named("end") = newend,
            Named("scores") = newscores,
            Named("cumulative.error") = cumulative_error);
}

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

List PTAProcessor::get_result() const {
    return List::create(
            Named("start") = start,
            Named("end") = end,
            Named("scores") = scores);
}

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
        stop("Intervals should be sorted and non-overlapping.");
    }
    return distance <= adjacency_treshold;
}

NumericVector PTAProcessor::merged_scores(int i, int j) const {
    NumericVector scores_i = const_cast<PTAProcessor*>(this)->scores(i, _);
    NumericVector scores_j = const_cast<PTAProcessor*>(this)->scores(j, _);
    return (length(i) * scores_i + length(j) * scores_j)
         / (length(i) + length(j));
}

double PTAProcessor::key(int heap, int nodeid) const {
    int previd = nodeid;
    for (int i = -1; i < heap; ++i) {
        if (previd == first_node) return INFINITY;
        previd = nodes[previd].prev;
    }

    switch (mode) {
        case PTA_MODE_NORMAL:
            return dsim(previd, nodeid);
        case PTA_MODE_CORRELATION:
            {
                double cor = correlation(previd, nodeid);
                return cor == 0 ? INFINITY : 1 / abs(cor);
            }
        default:
            stop("Bad PTA mode.");
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
        node.alive = false;
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

inline double sum_sq(const NumericVector x) {
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
        SEXP count_,
        SEXP error_,
        SEXP adjacency_treshold_,
        SEXP skip_,
        SEXP mode_,
        SEXP correlation_bound_)
    : start(start_), end(end_), scores(scores_)
{
    mode = as<int>(mode_);

    count_bound = as<int>(count_);
    error_bound = as<double>(error_);
    adjacency_treshold = as<double>(adjacency_treshold_);
    if (count_bound > 1) {
        maximum_error = INFINITY;
    } else if (mode == PTA_MODE_CORRELATION) {
        maximum_error = INFINITY;
        correlation_bound = as<double>(correlation_bound_);
    } else {
        // calculate maximum error
        NumericVector score1 = scores(0, _);
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
                    NumericVector diff = score1 - scores(j, _);
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
        node.alive = true;
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
        assert(is_heap(heaps[heap].begin(), heaps[heap].end(), greaters[heap]));
        for (int i = 0; i < node_count; ++i) {
            nodes[heaps[heap][i]].positions[heap] = i;
        }
    }
}

double PTAProcessor::dsim(int i, int j) const {
    if (!adjacent(i, j)) return INFINITY;

    const NumericVector z = merged_scores(i, j);
    NumericVector diff_i = z - const_cast<PTAProcessor *>(this)->scores(i, _);
    NumericVector diff_j = z - const_cast<PTAProcessor *>(this)->scores(j, _);
    return length(i) * sum_sq(diff_i) + length(j) * sum_sq(diff_j);
}

double PTAProcessor::correlation(int x, int y) const {
    if (!adjacent(x, y)) return 0;

    NumericVector scores_x = const_cast<PTAProcessor *>(this)->scores(x, _);
    NumericVector scores_y = const_cast<PTAProcessor *>(this)->scores(y, _);
    int n = scores_x.size();

    double mean_x = accumulate(scores_x.begin(), scores_x.end(), 0.0) / n;
    double mean_y = accumulate(scores_y.begin(), scores_y.end(), 0.0) / n;

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

void PTAProcessor::run() {
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
           if (1 / minkey <= correlation_bound) {
               break;
           }
           Rprintf("Correlation: %f\n", 1/minkey);
        } else {
            if (cumulative_error + minkey > abs_error_bound) {
                break;
            }
            cumulative_error += minkey;
        }

        if (!merge(minheap, minnode)) {
            break;
        }
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

    start = newstart;
    end = newend;
    scores = newscores;
}

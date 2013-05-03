#include <algorithm>
#include <cassert>
#include "MultiPTAProcessor.h"
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

double MultiPTAProcessor::key(int heap, int nodeid) const {
    int previd = nodeid;
    for (int i = -1; i < heap; ++i) {
        if (previd == first_node) return INFINITY;
        previd = nodes[previd].prev;
    }
    return dsim(previd, nodeid);
}

int MultiPTAProcessor::parent(int i) const {
    return (i - 1)/2;
}

int MultiPTAProcessor::left_child(int i) const {
    return 2 * i + 1;
}

int MultiPTAProcessor::right_child(int i) const {
    return 2 * i + 2;
}

void MultiPTAProcessor::heap_swap(int heap, int i, int j) {
    swap(heaps[heap][i], heaps[heap][j]);
    nodes[heaps[heap][i]].positions[heap] = i;
    nodes[heaps[heap][j]].positions[heap] = j;
}

void MultiPTAProcessor::heap_up(int heap, int i) {
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

void MultiPTAProcessor::heap_down(int heap, int i) {
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

void MultiPTAProcessor::heap_delete(int heap, int i) {
    heaps[heap][i] = heaps[heap].back();
    nodes[heaps[heap][i]].positions[heap] = i;
    heaps[heap].pop_back();

    if (i > 0 && greaters[heap](heaps[heap][parent(i)], heaps[heap][i])) {
        heap_up(heap, i);
    } else {
        heap_down(heap, i);
    }
}

void MultiPTAProcessor::heap_insert(int heap, int nodeid) {
    heaps[heap].push_back(nodeid);
    int i = heaps[heap].size() - 1;
    nodes[nodeid].positions[heap] = i;
    heap_up(heap, i);
}

void MultiPTAProcessor::update_node(int nodeid) {
    if (nodeid == -1) return;
    Node& node = nodes[nodeid];
    FOR_EACH_HEAP(heap) {
        heap_delete(heap, node.positions[heap]);
        node.keys[heap] = key(heap, node.id);
        heap_insert(heap, node.id);
    }
}

bool MultiPTAProcessor::merge(int minheap, int minnode) {
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

    score[repl.id] = merged_score(repl.id, minnode);
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

MultiPTAProcessor::MultiPTAProcessor(SEXP start_, SEXP end_, SEXP score_, SEXP count_, SEXP error_, SEXP adjacency_treshold_, SEXP skip_) :
    BasePTAProcessor(start_, end_, score_, count_, error_, adjacency_treshold_)
{
    nheaps = as<int>(skip_);
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

double MultiPTAProcessor::dsim(int i, int j) const {
    if (adjacent(i, j)) {
        const double z = merged_score(i, j);
        return length(i) * pow(z - score[i], 2) + length(j) * pow(z - score[j], 2);
    } else {
        return INFINITY;
    }
}

void MultiPTAProcessor::run() {
    double abs_error_bound = error_bound * maximum_error;
    double cumulative_error = 0;
    while (node_count > count_bound) {
        int minheap = 0;
        int minnode = heaps[0].front();
        double minerror = key(minheap, minnode);
        for (int i = 1; i < nheaps; ++i) {
            if (key(i, heaps[i].front()) < minerror) {
                minheap = i;
                minnode = heaps[i].front();
                minerror = key(minheap, minnode);
            }
        }

        if (cumulative_error + minerror > abs_error_bound) {
            break;
        }

        if (!merge(minheap, minnode)) {
            break;
        }

        cumulative_error += minerror;
    }

    NumericVector newstart(node_count);
    NumericVector newend(node_count);
    NumericVector newscore(node_count);
    int nodeid = first_node;
    int i = 0;
    while (nodeid >= 0) {
        newstart[i] = start[nodeid];
        newend[i] = end[nodeid];
        newscore[i] = score[nodeid];
        nodeid = nodes[nodeid].next;
        ++i;
    }

    start = newstart;
    end = newend;
    score = newscore;
}

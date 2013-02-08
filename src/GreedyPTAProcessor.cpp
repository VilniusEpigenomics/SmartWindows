#include <algorithm>
#include <cassert>
#include "GreedyPTAProcessor.h"
using namespace Rcpp;
using namespace std;

int GreedyPTAProcessor::peek() const {
    return heap.front();
}

double GreedyPTAProcessor::key(int nodeid) const {
    if (nodeid == first_node) {
        return INFINITY;
    } else {
        return dsim(nodes[nodeid].prev, nodeid);
    }
}

int GreedyPTAProcessor::parent(int i) {
    return (i - 1)/2;
}

int GreedyPTAProcessor::left_child(int i) {
    return 2 * i + 1;
}

int GreedyPTAProcessor::right_child(int i) {
    return 2 * i + 2;
}

void GreedyPTAProcessor::heap_swap(int i, int j) {
    swap(heap[i], heap[j]);
    nodes[heap[i]].position = i;
    nodes[heap[j]].position = j;
}

void GreedyPTAProcessor::heap_up(int i) {
    while (i > 0) {
        int par = parent(i);
        if (greater(heap[par], heap[i])) {
            heap_swap(i, par);
            i = par;
        } else {
            break;
        }
    }
}

void GreedyPTAProcessor::heap_down(int i) {
    while (true) {
        int smallest = i;
        int right = right_child(i);
        if (right < static_cast<int>(heap.size()) && greater(heap[smallest], heap[right])) {
            smallest = right;
        }
        int left = left_child(i);
        if (left < static_cast<int>(heap.size()) && greater(heap[smallest], heap[left])) {
            smallest = left;
        }

        if (smallest == i) {
            break;
        } else {
            heap_swap(i, smallest);
            i = smallest;
        }
    }
}

void GreedyPTAProcessor::heap_delete(int i) {
    heap[i] = heap.back();
    nodes[heap[i]].position = i;
    heap.pop_back();

    if (i > 0 && greater(heap[parent(i)], heap[i])) {
        heap_up(i);
    } else {
        heap_down(i);
    }
}

void GreedyPTAProcessor::heap_insert(int nodeid) {
    heap.push_back(nodeid);
    int i = heap.size() - 1;
    nodes[nodeid].position = i;
    heap_up(i);
}

void GreedyPTAProcessor::erase_node(int nodeid) {
    const Node& node = nodes[nodeid];

    if (nodeid == first_node) {
        first_node = node.next;
    } else {
        nodes[node.prev].next = node.next;
    }

    if (nodeid == last_node) {
        last_node = node.prev;
    } else {
        nodes[node.next].prev = node.prev;
    }

    --node_count;
}

bool GreedyPTAProcessor::merge() {
    assert(is_heap(heap.begin(), heap.end(), greater));
    if (heap.size() <= 1) return false;
    int topid = peek();
    const Node& top = nodes[topid];
    if (top.key == INFINITY) return false;

    score[top.prev] = merged_score(top.prev, topid);
    end[top.prev] = end[topid];
    heap_delete(0);
    assert(is_heap(heap.begin(), heap.end(), greater));
    erase_node(topid);

    if (top.next != -1) {
        Node& next = nodes[top.next];
        next.key = key(next.id);
        heap_delete(next.position);
        assert(is_heap(heap.begin(), heap.end(), greater));
        heap_insert(next.id);
        assert(is_heap(heap.begin(), heap.end(), greater));
    }


    return true;
}

GreedyPTAProcessor::GreedyPTAProcessor(SEXP start_, SEXP end_, SEXP score_, SEXP count_, SEXP error_) :
    BasePTAProcessor(start_, end_, score_, count_, error_),
    greater(*this)
{
    node_count = size();
    nodes.resize(node_count);
    heap.resize(node_count);
    first_node = 0;
    last_node = node_count - 1;
    for (int i = 0; i < node_count; ++i) {
        Node& node = nodes[i];
        node.prev = i - 1;
        node.id = i;
        node.next = i + 1;
        node.key = key(node.id);
        heap[i] = i;
    }
    nodes[last_node].next = -1;
    make_heap(heap.begin(), heap.end(), greater);
    assert(is_heap(heap.begin(), heap.end(), greater));
    for (int i = 0; i < node_count; ++i) {
        nodes[heap[i]].position = i;
    }
}

double GreedyPTAProcessor::dsim(int i, int j) const {
    if (adjacent(i, j)) {
        const double z = merged_score(i, j);
        return length(i) * pow(z - score[i], 2) + length(j) * pow(z - score[j], 2);
    } else {
        return INFINITY;
    }
}

void GreedyPTAProcessor::run() {
    while (node_count > count_bound) {
        if (!merge()) {
            break;
        }
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

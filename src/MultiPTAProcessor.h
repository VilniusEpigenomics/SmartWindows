#ifndef pta_MultiPTAProcessor_h
#define pta_MultiPTAProcessor_h
#include <Rcpp.h>
#include <vector>
#include <cassert>
#include "BasePTAProcessor.h"

class MultiPTAProcessor : public BasePTAProcessor {
    private:
        struct Node {
            int prev;
            int next;
            int id;
            bool alive;
            std::vector<int> positions;
            std::vector<double> keys;
        };

        struct NodeGreater {
            const MultiPTAProcessor* processor;
            int heap;
            NodeGreater() : processor(NULL), heap(-1) {}
            void set(const MultiPTAProcessor* p, int h) {
                processor = p;
                heap = h;
            }
            double operator()(int a, int b) {
                assert(processor);
                return processor->nodes[a].keys[heap] > processor->nodes[b].keys[heap];
            }
        };

        std::vector<NodeGreater> greaters;
        int nheaps;
        std::vector< std::vector<int> > heaps;
        std::vector<Node> nodes;
        int node_count;
        int first_node;
        int last_node;

        double dsim(int i, int j) const;

        int parent(int i) const;
        int left_child(int i) const;
        int right_child(int i) const;

        double key(int heap, int nodeid) const;
        void heap_swap(int heap, int i, int j);
        void heap_up(int heap, int i);
        void heap_down(int heap, int i);
        void heap_delete(int heap, int i);
        void heap_insert(int heap, int nodeid);

        void update_node(int nodeid);
        bool merge(int heap, int node);

    public:
        MultiPTAProcessor(SEXP start_, SEXP end_, SEXP score_, SEXP count_, SEXP error_, SEXP adjacency_treshold_, SEXP skip_);
        void run();
};

#endif

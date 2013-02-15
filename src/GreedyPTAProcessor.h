#ifndef pta_GreedyPTAProcessor_h
#define pta_GreedyPTAProcessor_h
#include <Rcpp.h>
#include <vector>
#include "BasePTAProcessor.h"

class GreedyPTAProcessor : public BasePTAProcessor {
    private:
        struct Node {
            int prev;
            int next;
            int id;
            int position;
            double key;
        };

        struct NodeGreater {
            const GreedyPTAProcessor& processor;
            NodeGreater(const GreedyPTAProcessor& p) : processor(p) {}
            double operator()(int a, int b) {
                return processor.nodes[a].key > processor.nodes[b].key;
            }
        };

        NodeGreater greater;
        std::vector<int> heap;
        std::vector<Node> nodes;
        int node_count;
        int first_node;
        int last_node;

        double dsim(int i, int j) const;
        double key(int nodeid) const;

        int parent(int i);
        int left_child(int i);
        int right_child(int i);
        void heap_swap(int i, int j);
        void heap_up(int i);
        void heap_down(int i);
        void heap_delete(int i);
        void heap_insert(int nodeid);

        int peek() const;
        void erase_node(int nodeid);
        bool merge();

    public:
        GreedyPTAProcessor(SEXP start_, SEXP end_, SEXP score_, SEXP count_, SEXP error_, SEXP adjacency_treshold_);
        void run();
};

#endif

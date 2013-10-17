#ifndef PTA_PTAProcessor_h
#define PTA_PTAProcessor_h
#include <Rcpp.h>
#include <vector>
#include <cassert>

#define PTA_MODE_NORMAL 0
#define PTA_MODE_CORRELATION 1

class __attribute__((visibility("hidden"))) PTAProcessor {
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
            const PTAProcessor* processor;
            int heap;
            NodeGreater() : processor(NULL), heap(-1) {}
            void set(const PTAProcessor* p, int h) {
                processor = p;
                heap = h;
            }
            double operator()(int a, int b) {
                assert(processor);
                return processor->nodes[a].keys[heap] > processor->nodes[b].keys[heap];
            }
        };

        const Rcpp::NumericVector original_start;
        const Rcpp::NumericVector original_end;
        const Rcpp::NumericMatrix original_scores;

        Rcpp::NumericVector start;
        Rcpp::NumericVector end;
        Rcpp::NumericMatrix scores;

        int count_bound;
        double error_bound;
        double adjacency_threshold;
        double correlation_bound;
        bool correlation_spearman;
        bool correlation_absolute;

        int minimum_count;
        double maximum_error;

        int mode;

        inline int size() const { return start.size(); }
        double length(int interval) const;
        bool adjacent(int i, int j) const;
        Rcpp::NumericVector merged_scores(int i, int j) const;

        std::vector<NodeGreater> greaters;
        int nheaps;
        std::vector< std::vector<int> > heaps;
        std::vector<Node> nodes;
        int node_count;
        int first_node;
        int last_node;

        double dsim(int i, int j) const;
        double correlation(int i, int j) const;

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
        PTAProcessor(const Rcpp::List arguments);
        Rcpp::List run();
};

#endif

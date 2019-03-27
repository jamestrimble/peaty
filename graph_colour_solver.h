#ifndef GRAPH_COLOUR_SOLVER_H
#define GRAPH_COLOUR_SOLVER_H

#include <atomic>
#include <vector>

struct ColouringGraph {
    int n;
    std::vector<std::vector<bool>> adj_matrix;
    std::vector<std::vector<int>> adjlist;

    ColouringGraph(int n) : n(n), adj_matrix(n, std::vector<bool>(n)), adjlist(n) {
    }

    void add_edge(int v, int w)
    {
        adj_matrix[v][w] = true;
        adj_matrix[w][v] = true;
    }

    struct ColouringGraph induced_subgraph(std::vector<int> & vv) const;

    void make_adjacency_lists();
};


struct Solution {
    int size = 0;
    std::vector<std::vector<int>> vtx_colour;
    Solution(int n, int f) : vtx_colour(n, std::vector<int>(f)) {}
};


//void solve(struct ColouringGraph *g, unsigned long long *expand_call_count, unsigned long long expand_call_limit,
//        bool quiet, struct Solution *incumbent, int num_colours, int f);
//
//bool is_solution_valid(struct ColouringGraph *g, struct Solution *solution,
//        int num_colours);

int find_colouring_number(const ColouringGraph & g, int f, std::atomic_bool & terminate_early);

class ColouringNumberFinder
{
    const ColouringGraph & g;
    ColouringGraph sorted_graph;
    int f;
    int current_target_num_colours = 0;
    unsigned rng_seed = 0;
    unsigned long long search_node_count = 0;
    unsigned long long local_search_node_limit = 1000;
    int colouring_number = -1;

public:
    unsigned long long get_search_node_count();
    int get_colouring_number();
    void search();
    ColouringNumberFinder(const ColouringGraph & g, int f);
};

#endif

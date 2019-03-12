#include <vector>

struct ColouringGraph {
    int n;
    std::vector<std::vector<bool>> adj_matrix;
    std::vector<std::vector<int>> adjlist;

    ColouringGraph(int n) : n(n), adj_matrix(n, std::vector<bool>(n)), adjlist(n) {
    }
};

void add_edge(struct ColouringGraph & g, int v, int w);

struct ColouringGraph induced_subgraph(const ColouringGraph *g, std::vector<int> & vv);

void make_adjacency_lists(struct ColouringGraph *g);


struct Solution {
    int size;
    int capacity;
    int *vtx_colour;
};

void init_Solution(struct Solution *l, int capacity);

void destroy_Solution(struct Solution *l);


//void solve(struct ColouringGraph *g, unsigned long long *expand_call_count, unsigned long long expand_call_limit,
//        bool quiet, struct Solution *incumbent, int num_colours, int f);
//
//bool is_solution_valid(struct ColouringGraph *g, struct Solution *solution,
//        int num_colours);

int find_colouring_number(const ColouringGraph & g, int f);

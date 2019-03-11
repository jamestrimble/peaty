#include <vector>

struct Graph {
    int n;
    std::vector<std::vector<bool>> adj_matrix;
    std::vector<std::vector<int>> adjlist;

    Graph(int n) : n(n), adj_matrix(n, std::vector<bool>(n)), adjlist(n) {
    }
};

void add_edge(struct Graph & g, int v, int w);

struct Graph induced_subgraph(const Graph *g, std::vector<int> & vv);

void make_adjacency_lists(struct Graph *g);


struct Solution {
    int size;
    int capacity;
    int *vtx_colour;
};

void init_Solution(struct Solution *l, int capacity);

void destroy_Solution(struct Solution *l);


//void solve(struct Graph *original_g, unsigned long long *expand_call_count, unsigned long long expand_call_limit,
//        bool quiet, struct Solution *incumbent, int num_colours, int f);
//
//bool is_solution_valid(struct Graph *original_g, struct Solution *solution,
//        int num_colours);

int find_colouring_number(const Graph & g, int f);

struct Graph {
    int n;
    int *degree;
    bool **adj_matrix;
    int *adjlist_len;
    int **adjlist;
};

void remove_edge(struct Graph *g, int v, int w);

void add_edge(struct Graph *g, int v, int w);

struct Graph *induced_subgraph(struct Graph *g, int *vv, int vv_len);

struct Graph *new_graph(int n);

void free_graph(struct Graph *g);

void make_adjacency_lists(struct Graph *g);

void free_adjacency_lists(struct Graph *g);



struct Solution {
    int size;
    int *vtx_colour;
};

void init_Solution(struct Solution *l, int capacity);

void destroy_Solution(struct Solution *l);



void solve(struct Graph *original_g, long *expand_call_count,
        bool quiet, struct Solution *incumbent, int num_colours);

bool is_solution_valid(struct Graph *original_g, struct Solution *solution,
        int num_colours);


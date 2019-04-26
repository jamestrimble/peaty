#define _POSIX_SOURCE

#include <argp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <functional>
#include <mutex>

#include "graph.h"
#include "sparse_graph.h"
#include "util.h"
#include "sequential_solver.h"
#include "params.h"

using std::atomic;
using std::condition_variable;
using std::cv_status;
using std::function;
using std::mutex;
using std::unique_lock;

using std::chrono::seconds;
using std::chrono::steady_clock;
using std::chrono::system_clock;
using std::chrono::duration_cast;
using std::chrono::milliseconds;

static char doc[] = "Find a maximum clique in a graph in DIMACS format";
static char args_doc[] = "";
static struct argp_option options[] = {
    {"quiet", 'q', 0, 0, "Quiet output"},
    {"unweighted-sort", 'u', 0, 0, "Unweighted ordering (only applies to certain algorithms)"},
    {"colouring-variant", 'c', "VARIANT", 0, "For algorithms 0 and 5, which type of colouring?"},
    {"algorithm", 'a', "NUMBER", 0, "Algorithm number"},
    {"max-sat-level", 'm', "LEVEL", 0, "Level of MAXSAT reasoning; default=2"},
    {"num-threads", 't', "NUMBER", 0, "Number of threads (for parallel algorithms only)"},
    {"file-format", 'f', "FORMAT", 0, "File format (DIMACS, MTX or EDGES)"},
    { 0 }
};

enum class FileFormat
{
    Dimacs,
    Pace
};

static struct {
    bool quiet;
    bool unweighted_sort;
    int colouring_variant = 3;
    int algorithm_num;
    int max_sat_level = -1;
    int num_threads = 1;
    FileFormat file_format = FileFormat::Pace;
} arguments;

static error_t parse_opt (int key, char *arg, struct argp_state *state)
{
    switch (key) {
        case 'q':
            arguments.quiet = true;
            break;
        case 'u':
            arguments.unweighted_sort = true;
            break;
        case 'c':
            arguments.colouring_variant = atoi(arg);
            break;
        case 'a':
            arguments.algorithm_num = atoi(arg);
            break;
        case 'm':
            arguments.max_sat_level = atoi(arg);
            break;
        case 't':
            arguments.num_threads = atoi(arg);
            break;
        case 'f':
            if (!strcmp(arg, "PACE") || !strcmp(arg, "pace"))
                arguments.file_format = FileFormat::Pace;
            else if (!strcmp(arg, "DIMACS") || !strcmp(arg, "dimacs"))
                arguments.file_format = FileFormat::Dimacs;
            break;
        case ARGP_KEY_ARG:
//            argp_usage(state);
            break;
        case ARGP_KEY_END:
            break;
        default: return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

static struct argp argp = { options, parse_opt, args_doc, doc };

/******************************************************************************/

// Checks if a set of vertices induces a clique
bool check_vertex_cover(const SparseGraph & g, const vector<int> & vc) {
    // TODO: make sure the graph passed in here isn't modified from the original
    vector<bool> in_vc(g.n);
    for (int v : vc) {
        in_vc[v] = true;
    }
    for (unsigned i=0; i<g.n; i++) {
        if (g.vertex_has_loop[i] && !in_vc[i]) {
            std::cerr << "Vertex " << i << " has a loop but is not in the vertex cover!" << std::endl;
            return false;
        }
    }
    for (unsigned i=0; i<g.n; i++) {
        if (!in_vc[i]) {
            for (unsigned v : g.adjlist[i]) {
                if (!in_vc[v]) {
                    std::cerr << "Edge " << i << "," << v << " is uncovered!" << std::endl;
                    return false;
                }
            }
        }
    }
    return true;
}

struct Result
{
    VtxList vertex_cover;
    long search_node_count;
    Result(const SparseGraph & g) : vertex_cover(g.n), search_node_count(0) {}
};

struct Reduction
{
    virtual void unwind(vector<bool> & in_cover) = 0;
    virtual ~Reduction() {}
};

// a list of tuples (v, w, x), where v is the vertex of deg 2,
// w is the kept neighbour and x is the removed neighbour
struct Deg2Reduction : public Reduction
{
    int v;
    int w;
    int x;

    Deg2Reduction(int v, int w, int x) : v(v), w(w), x(x) {}

    void unwind(vector<bool> & in_cover)
    {
        if (in_cover[w]) {
            in_cover[x] = true;
        } else {
            in_cover[v] = true;
        }
    }

    ~Deg2Reduction() {}
};

// a list of tuples (v, w, x, y), where v is the vertex of deg 3,
// and w and x are the kept neighbours
struct FunnelReduction : public Reduction
{
    int v;
    vector<int> ww;
    int y;

    FunnelReduction(int v, vector<int> ww, int y) : v(v), ww(ww), y(y) {}

    bool all_ww_are_in_cover(vector<bool> & in_cover)
    {
        for (int w : ww)
            if (!in_cover[w])
                return false;
        return true;
    }

    void unwind(vector<bool> & in_cover)
    {
        if (all_ww_are_in_cover(in_cover)) {
            in_cover[y] = true;
        } else {
            in_cover[v] = true;
        }
    }

    ~FunnelReduction() {}
};

struct BowTieReduction : public Reduction
{
    int v;
    int a;
    int b;
    int c;
    int d;

    BowTieReduction(int v, int a, int b, int c, int d) : v(v), a(a), b(b), c(c), d(d) {}

    void unwind(vector<bool> & in_cover)
    {
        if (!in_cover[a]) {
            in_cover[c] = false;
            in_cover[v] = true;
        } else if (!in_cover[b]) {
            in_cover[d] = false;
            in_cover[v] = true;
        } else if (!in_cover[c]) {
            in_cover[b] = false;
            in_cover[v] = true;
        } else if (!in_cover[d]) {
            in_cover[a] = false;
            in_cover[v] = true;
        } else {
            in_cover[v] = false;
        }
    }

    ~BowTieReduction() {}
};

bool isolated_vertex_removal(SparseGraph & g,
        vector<bool> & in_cover, vector<bool> & deleted)
{
    bool made_a_change = false;
    vector<int> neighbours;
    for (unsigned v=0; v<g.n; v++) {
        if (!deleted[v]) {
            neighbours.clear();
            for (int w : g.adjlist[v])
                if (!deleted[w])
                    neighbours.push_back(w);
            if (g.vv_are_clique(neighbours)) {
                deleted[v] = true;
                made_a_change = true;
                for (int w : neighbours) {
                    deleted[w] = true;
                    in_cover[w] = true;
                    for (int u : g.adjlist[w]) {
                        auto & lst = g.adjlist[u];
                        lst.erase(std::find(lst.begin(), lst.end(), w));
                    }
                    g.adjlist[w].clear();
                }
            }
        }
    }
    return made_a_change;
}

// Remove w from the adjacency list of v.
// It is the caller's responsibility to ensure that v gets removed
// from the adjacency list of w.
auto remove_from_adj_list(SparseGraph & g, int v, int w)
{
    auto & lst = g.adjlist[v];
    lst.erase(std::find(lst.begin(), lst.end(), w));
}

bool vertex_folding(SparseGraph & g, vector<bool> & deleted,
        vector<std::unique_ptr<Reduction>> & reductions)
{
    bool made_a_change = false;

    for (unsigned v=0; v<g.n; v++) {
        if (g.adjlist[v].size() == 2) {
            int w = g.adjlist[v][0];
            int x = g.adjlist[v][1];

            // for this reduction, w and x must not be adjacent
            if (!g.has_edge(x, w)) {
                remove_from_adj_list(g, w, v);
                remove_from_adj_list(g, x, v);
                for (int u : g.adjlist[x]) {
                    remove_from_adj_list(g, u, x);
                    if (!g.has_edge(u, w))
                        g.add_edge(w, u);
                }
                g.adjlist[v].clear();
                g.adjlist[x].clear();
                deleted[v] = true;
                deleted[x] = true;
                reductions.push_back(std::make_unique<Deg2Reduction>(v, w, x));
                made_a_change = true;
            }
        }
    }
    return made_a_change;
}

bool has_any_edge(SparseGraph & g, int v, vector<int> & ww)
{
    for (int w : ww)
        if (w != v && g.has_edge(v, w))
                return true;
    return false;
}

int num_edges(SparseGraph & g, vector<int> & vv)
{
    int retval = 0;
    for (int v : vv) {
        for (int w : vv) {
            if (w > v && g.has_edge(v, w)) {
                ++retval;
            }
        }
    }
    return retval;
}

// TODO: maybe do the normal, more general version of funnel
// in which y can be adjacent so some of ww
bool do_funnel_reductions(SparseGraph & g, vector<bool> & in_cover, vector<bool> & deleted,
        vector<std::unique_ptr<Reduction>> & reductions)
{
    bool made_a_change = false;

    // Try to find vertex v with neighbours ww and y, such that ws are a clique
    // and y is not adjacent to any member of ww.
    // We can then delete v and y, adding the neighbours
    // of y other than x to the adjacency list of each vertex in ww.
    for (unsigned v=0; v<g.n; v++) {
        int lst_sz = g.adjlist[v].size();
        if (lst_sz >= 3) {
            if (num_edges(g, g.adjlist[v]) != (lst_sz-1) * (lst_sz-2) / 2)
                continue;

            for (int y : g.adjlist[v]) {
                if (!has_any_edge(g, y, g.adjlist[v])) {
                    vector<int> ww;
                    for (int w : g.adjlist[v])
                        if (w != y)
                            ww.push_back(w);

                    for (int w : ww) {
                        remove_from_adj_list(g, w, v);
                    }
                    remove_from_adj_list(g, y, v);
                    for (int u : g.adjlist[y]) {
                        remove_from_adj_list(g, u, y);
                        for (int w : ww)
                            if (!g.has_edge(u, w))
                                g.add_edge(u, w);
                    }
                    g.adjlist[v].clear();
                    g.adjlist[y].clear();
                    deleted[v] = true;
                    deleted[y] = true;
                    reductions.push_back(std::make_unique<FunnelReduction>(v, ww, y));
                    made_a_change = true;

                    break;
                }
            }
        }
    }
    return made_a_change;
}


void check_adj_list_integrity(SparseGraph & g)
{
    for (unsigned v=0; v<g.n; v++) {
        std::vector<bool> a(g.n);
        for (int w : g.adjlist[v]) {
            if (a[w]) {
                std::cout << "Duplicate edge" << std::endl;
                exit(1);
            }
            a[w] = true;
            auto it = std::find(g.adjlist[w].begin(), g.adjlist[w].end(), v);
            if (it == g.adjlist[w].end()) {
                std::cout << "Graph error" << std::endl;
                exit(1);
            }
        }
    }
}


// If there is a vertex v with a neighbour x who is adjacent to all of v's other neighbours,
// it's safe to assume that x is in the vertex cover.
bool do_domination_reductions(SparseGraph & g, vector<bool> & in_cover, vector<bool> & deleted,
        vector<std::unique_ptr<Reduction>> & reductions)
{
    bool made_a_change = false;

    for (unsigned v=0; v<g.n; v++) {
        if (g.adjlist[v].size() > 2) {
            for (int w : g.adjlist[v]) {
                int num_edges = 0;
                for (int x : g.adjlist[v]) {
                    if (x != w) {
                        if (g.has_edge(x, w)) {
                            ++num_edges;
                        } else {
                            // for efficiency.  Maybe this can be tidied?
                            break;
                        }
                    }
                }
                if (num_edges == g.adjlist[v].size() - 1) {
                    for (int u : g.adjlist[w]) {
                        remove_from_adj_list(g, u, w);
                    }
                    g.adjlist[w].clear();
                    deleted[w] = true;
                    in_cover[w] = true;
                    made_a_change = true;
                    break;
                }
            }
        }
    }
    return made_a_change;
}

bool are_bow_tie(SparseGraph & g, const vector<int> & adjlist)
{
    for (int v : adjlist) {
        int num_edges = 0;
        for (int w : adjlist) {
            if (v != w) {
                if (g.has_edge(v, w)) {
                    ++num_edges;
                    if (num_edges > 1) {
                        break;
                    }
                }
            }
        }
        if (num_edges != 1) {
            return false;
        }
    }
    return true;
}

bool do_bow_tie_reductions(SparseGraph & g, vector<bool> & in_cover, vector<bool> & deleted,
        vector<std::unique_ptr<Reduction>> & reductions)
{
    bool made_a_change = false;

    for (unsigned v=0; v<g.n; v++) {
        if (g.adjlist[v].size() == 4) {
            auto & lst_v = g.adjlist[v];
            if (are_bow_tie(g, lst_v)) {
                int a = lst_v[0];
                int b = lst_v[1];
                int c = lst_v[2];
                int d = lst_v[3];
                if (g.has_edge(a, c)) {
                    std::swap(b, c);
                } else if (g.has_edge(a, d)) {
                    std::swap(b, d);
                }

                for (int u : g.adjlist[v]) {
                    remove_from_adj_list(g, u, v);
                }
                g.adjlist[v].clear();
                deleted[v] = true;

                vector<int> adjlist_a = g.adjlist[a];
                vector<int> adjlist_b = g.adjlist[b];
                vector<int> adjlist_c = g.adjlist[c];
                vector<int> adjlist_d = g.adjlist[d];
                for (int u : adjlist_c)
                    if (!g.has_edge(a, u))
                        g.add_edge(a, u);
                for (int u : adjlist_d)
                    if (!g.has_edge(b, u))
                        g.add_edge(b, u);
                for (int u : adjlist_b)
                    if (!g.has_edge(c, u))
                        g.add_edge(c, u);
                for (int u : adjlist_a)
                    if (!g.has_edge(d, u))
                        g.add_edge(d, u);

                reductions.push_back(std::make_unique<BowTieReduction>(v, a, b, c, d));
                made_a_change = true;
            }
        }
    }

    return made_a_change;
}

auto make_list_of_components(const SparseGraph & g) -> vector<vector<int>>
{
    vector<vector<int>> components;
    vector<bool> vertex_used(g.n);
    for (unsigned i=0; i<g.n; i++)
        if (g.adjlist[i].empty())
            vertex_used[i] = true;

    for (unsigned i=0; i<g.n; i++) {
        if (!vertex_used[i]) {
            components.push_back({int(i)});
            auto & component = components.back();
            vector<int> to_explore = {int(i)};
            vertex_used[i] = true;
            while (!to_explore.empty()) {
                int v = to_explore.back();
                to_explore.pop_back();
                for (int w : g.adjlist[v]) {
                    if (!vertex_used[w]) {
                        component.push_back(w);
                        to_explore.push_back(w);
                        vertex_used[w] = true;
                    }
                }
            }
        }
    }

    return components;
}

auto find_vertex_cover_of_subgraph(const SparseGraph & g, vector<int> component,
        const Params & params) -> vector<int>
{
    std::sort(component.begin(), component.end());
    SparseGraph subgraph = g.induced_subgraph<SparseGraph>(component);

//    for (unsigned v=0; v<subgraph.n; v++) {
//        for (int w : subgraph.adjlist[v]) {
//            for (int u : subgraph.adjlist[v]) {
//                std::cout << (subgraph.has_edge(w, u) ? "X " : ". ");
//            }
//            std::cout << std::endl;
//        }
//        std::cout << std::endl;
//    }
//    subgraph.print_dimacs_format();

    VtxList independent_set(g.n);
    long search_node_count = 0;

    sequential_mwc(subgraph, params, independent_set, search_node_count);

    vector<bool> vtx_is_in_ind_set(subgraph.n);
    for (int v : independent_set.vv) {
        vtx_is_in_ind_set[v] = true;
    }
    vector<int> vertex_cover;
    for (unsigned i=0; i<component.size(); i++) {
        if (!vtx_is_in_ind_set[i]) {
            vertex_cover.push_back(component[i]);
        }
    }
    return vertex_cover;
}

auto mwc(SparseGraph g, const Params & params) -> Result
{
    vector<bool> deleted = g.vertex_has_loop;
    vector<bool> in_cover = g.vertex_has_loop;
    g.remove_edges_incident_to_loopy_vertices();

    vector<std::unique_ptr<Reduction>> reductions;

    while (true) {
        bool a = isolated_vertex_removal(g, in_cover, deleted);
        bool b = do_domination_reductions(g, in_cover, deleted, reductions);
        bool c = vertex_folding(g, deleted, reductions);
        bool d = do_funnel_reductions(g, in_cover, deleted, reductions);
        bool e = false;  //do_bow_tie_reductions(g, in_cover, deleted, reductions);
        if (!a && !b && !c && !d && !e)
            break;
    };
    check_adj_list_integrity(g);

    vector<vector<int>> components = make_list_of_components(g);

//    for (auto & component : components) {
//        std::cout << "A_COMPONENT";
//        for (int v : component) {
//            std::cout << " " << v;
//        }
//        std::cout << std::endl;
//    }
//    std::cout << "END_COMPONENTS" << std::endl;

    Result result(g);
    for (auto & component : components) {
        std::cout << "c COMPONENT " << component.size() << std::endl;
        auto vertex_cover_of_subgraph = find_vertex_cover_of_subgraph(g, component, params);
        for (int v : vertex_cover_of_subgraph) {
            in_cover[v] = true;
        }
    }

    while (!reductions.empty()) {
        reductions.back()->unwind(in_cover);
        reductions.pop_back();
    }

    for (unsigned v=0; v<g.n; v++) {
        if (in_cover[v]) {
            result.vertex_cover.push_vtx(v, 1);
        }
    }

    return result;
}

int main(int argc, char** argv) {
    argp_parse(&argp, argc, argv, 0, 0, 0);

    if (arguments.algorithm_num != 5)
        arguments.num_threads = 1;

    SparseGraph g =
            arguments.file_format==FileFormat::Pace ? readSparseGraphPaceFormat() :
                                                      readSparseGraph();

    Params params {arguments.colouring_variant, arguments.max_sat_level, arguments.algorithm_num,
            arguments.num_threads, arguments.quiet, arguments.unweighted_sort};

    Result result = mwc(g, params);

    // sort vertices in clique by index
    std::sort(result.vertex_cover.vv.begin(), result.vertex_cover.vv.end());

    std::cout << "s vc " << g.n << " " << result.vertex_cover.vv.size() << std::endl;
    for (int v : result.vertex_cover.vv)
        std::cout << (v+1) << std::endl;

//    printf("Stats: status program algorithm_number max_sat_level num_threads size weight nodes\n");
//    std::cout <<
//            (aborted ? "TIMEOUT" : "COMPLETED") << " " <<
//            argv[0] << " " <<
//            arguments.algorithm_num << " " <<
//            arguments.max_sat_level << " " <<
//            arguments.num_threads << " " <<
//            result.vertex_cover.vv.size() << " " <<
//            result.vertex_cover.total_wt <<  " " <<
//            result.search_node_count << std::endl;

    if (!check_vertex_cover(g, result.vertex_cover.vv))
        fail("*** Error: invalid solution\n");
}

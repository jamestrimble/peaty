#define _POSIX_SOURCE

#include "graph.h"
#include "sparse_graph.h"
#include "util.h"
#include "sequential_solver.h"
#include "params.h"

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
#include <thread>

using std::atomic;
using std::condition_variable;
using std::cv_status;
using std::function;
using std::mutex;
using std::thread;
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
    {"ind-set-upper-bound", 'i', "BOUND", 0, "An upper bound on the size of the max independent set"},
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
    int ind_set_upper_bound;
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
        case 'i':
            arguments.ind_set_upper_bound = atoi(arg);
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

// a list of tuples (v, w, x), where v is the vertex of deg 2,
// w is the kept neighbour and x is the removed neighbour
struct Deg2Reduction
{
    int v;
    int w;
    int x;
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

bool vertex_folding(SparseGraph & g,
        vector<bool> & in_cover, vector<bool> & deleted,
        vector<Deg2Reduction> & deg_2_reductions)
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
                deg_2_reductions.push_back({int(v), w, x});
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

auto mwc(SparseGraph & g, const Params & params) -> Result
{
    vector<bool> deleted = g.vertex_has_loop;
    vector<bool> in_cover = g.vertex_has_loop;
    g.remove_edges_incident_to_loopy_vertices();

    vector<Deg2Reduction> deg_2_reductions;

    while (true) {
        bool a = isolated_vertex_removal(g, in_cover, deleted);
        bool b = vertex_folding(g, in_cover, deleted, deg_2_reductions);
        if (!a && !b)
            break;
    };

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
        auto vertex_cover_of_subgraph = find_vertex_cover_of_subgraph(g, component, params);
        for (int v : vertex_cover_of_subgraph) {
            in_cover[v] = true;
        }
    }

    while (!deg_2_reductions.empty()) {
        auto r = deg_2_reductions.back();
        deg_2_reductions.pop_back();
        if (in_cover[r.w]) {
            in_cover[r.x] = true;
        } else {
            in_cover[r.v] = true;
        }
    }
//    sequential_mwc(g, params, result.clq, result.search_node_count);
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
            arguments.file_format==FileFormat::Pace ? fastReadSparseGraphPaceFormat() :
                                                      readSparseGraph();

    Params params {arguments.colouring_variant, arguments.max_sat_level, arguments.algorithm_num,
            arguments.num_threads, arguments.quiet, arguments.unweighted_sort, arguments.ind_set_upper_bound};

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

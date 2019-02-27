#include "bitset.h"
#include "root_node_processing.h"
#include "graph.h"
#include "util.h"
#include "params.h"

#include <algorithm>
#include <vector>

using std::vector;
class SimpleColourer {
    Graph & g;
    const Params params;

    vector<unsigned long long> to_colour;
    vector<unsigned long long> candidates;
    vector<long> residual_wt;
    vector<int> col_class;

public:
    SimpleColourer(Graph & g, const Params params)
            : g(g), params(params), to_colour(g.numwords),
              candidates(g.numwords), residual_wt(g.n)
    {
        col_class.reserve(g.n);
    }

    auto simple_colouring_bound(vector<unsigned long long> & P_bitset,
            vector<long> & bounds, vector<int> & P) -> long
    {
        int numwords = calc_numwords(P_bitset, g.numwords);
        if (numwords == 0)
            return false;

        copy_bitset(P_bitset, to_colour, numwords);
        residual_wt = g.weight;

        long bound = 0;
        int v;
        while ((v=first_set_bit(to_colour, numwords))!=-1) {
            long class_min_wt = residual_wt[v];
            col_class.clear();
            col_class.push_back(v);
            bitset_intersection(to_colour, g.bit_complement_nd[v], candidates, numwords);
            while ((v=first_set_bit(candidates, numwords))!=-1) {
                if (residual_wt[v] < class_min_wt)
                    class_min_wt = residual_wt[v];
                col_class.push_back(v);
                bitset_intersect_with(candidates, g.bit_complement_nd[v], numwords);
            }
            bound += class_min_wt;
            for (int w : col_class) {
                residual_wt[w] -= class_min_wt;
                unset_bit_if(to_colour, w, residual_wt[w] == 0);
                if(residual_wt[w] == 0) {
                    bounds.push_back(bound);
                    P.push_back(w);
                }
            }
        }
#ifdef PRINT_BOUND_AND_EXIT
        printf("%ld\n", bound);
        exit(0);
#endif
        return bound;
    }
};

class SimpleMWC {
    Graph & g;
    const Params params;
    SimpleColourer colourer;
    VtxList & incumbent;
    long search_node_count = 0;

    void expand(VtxList& C, vector<unsigned long long> & P_bitset, long & search_node_count)
    {
        ++search_node_count;
        if (params.abort->load())
            return;

        if (bitset_empty(P_bitset, g.numwords) && C.total_wt > incumbent.total_wt) {
            incumbent = C;
            if (!params.quiet)
                printf("New incumbent of weight %ld\n", incumbent.total_wt);
        }

        vector<long> bounds;
        vector<int> P;
        bounds.reserve(g.n);
        P.reserve(g.n);

        if (colourer.simple_colouring_bound(P_bitset, bounds, P) + C.total_wt > incumbent.total_wt) {
            vector<unsigned long long> new_P_bitset(g.numwords);

            switch (params.algorithm_num)
            {
            case 1:
                {
                    while (!P.empty() && C.total_wt + bounds.back() > incumbent.total_wt) {
                        int v = P.back();
                        unset_bit(P_bitset, v);
                        bitset_intersection_with_complement(P_bitset, g.bit_complement_nd[v], new_P_bitset, g.numwords);
                        C.push_vtx(v, g);
                        expand(C, new_P_bitset, search_node_count);
                        C.pop_vtx(g);
                        P.pop_back();
                        bounds.pop_back();
                    }
                }
                break;
            case 2:
                {
                    unsigned start_idx = P.size();
                    while (start_idx != 0 && C.total_wt + bounds[start_idx - 1] > incumbent.total_wt) {
                        --start_idx;
                        unset_bit(P_bitset, P[start_idx]);
                    }
                    for (unsigned idx = start_idx; idx < P.size(); idx++) {
                        int v = P[idx];
                        bitset_intersection_with_complement(P_bitset, g.bit_complement_nd[v], new_P_bitset, g.numwords);
                        C.push_vtx(v, g);
                        expand(C, new_P_bitset, search_node_count);
                        C.pop_vtx(g);
                        set_bit(P_bitset, v);
                    }
                }
                break;
            case 3:
                {
                    vector<unsigned long long> branch_vv_bitset(g.numwords);
                    unsigned start_idx = P.size();
                    while (start_idx != 0 && C.total_wt + bounds[start_idx - 1] > incumbent.total_wt) {
                        --start_idx;
                        set_bit(branch_vv_bitset, P[start_idx]);
                    }
                    bitset_intersect_with_complement(P_bitset, branch_vv_bitset, g.numwords);
                    int v;
                    while ((v=first_set_bit(branch_vv_bitset, g.numwords))!=-1) {
                        unset_bit(branch_vv_bitset, v);
                        bitset_intersection_with_complement(P_bitset, g.bit_complement_nd[v], new_P_bitset, g.numwords);
                        C.push_vtx(v, g);
                        expand(C, new_P_bitset, search_node_count);
                        set_bit(P_bitset, v);
                        C.pop_vtx(g);
                    }
                }
                break;
            case 4:
                {
                    vector<unsigned long long> branch_vv_bitset(g.numwords);
                    unsigned start_idx = P.size();
                    while (start_idx != 0 && C.total_wt + bounds[start_idx - 1] > incumbent.total_wt) {
                        --start_idx;
                        set_bit(branch_vv_bitset, P[start_idx]);
                    }
                    int v;
                    while ((v=last_set_bit(branch_vv_bitset, g.numwords))!=-1) {
                        unset_bit(branch_vv_bitset, v);
                        unset_bit(P_bitset, v);
                        bitset_intersection_with_complement(P_bitset, g.bit_complement_nd[v], new_P_bitset, g.numwords);
                        C.push_vtx(v, g);
                        expand(C, new_P_bitset, search_node_count);
                        C.pop_vtx(g);
                    }
                }
                break;
            }
        }
    }

public:
    SimpleMWC(Graph & g, const Params params, VtxList & incumbent)
            : g(g), params(params), colourer(g, params), incumbent(incumbent)
    {
    }

    auto run(VtxList & C, long & search_node_count) -> void
    {
        vector<unsigned long long> P_bitset(g.numwords, 0);
        set_first_n_bits(P_bitset, g.n);
        expand(C, P_bitset, search_node_count);
    }
};

vector<int> calc_degs(const SparseGraph & g) {
    vector<int> degs(g.n);
    for (int v=0; v<int(g.n); v++)
        for (auto w : g.adjlist[v])
            ++degs[w];

    return degs;
}

auto reduce_and_reverse_unweighted_cp_order(vector<int> & vv, const SparseGraph & g) -> void
{
    if (vv.empty())
        return;

    auto residual_degs = calc_degs(g);

    unsigned new_sz = vv.size();

    for (unsigned i=vv.size(); i--; ) {
        // find vertex with lowest residual_deg
        unsigned best_v_pos = i;
        int best_deg = residual_degs[vv[i]];
        for (unsigned j=i; j--; ) {
            int v = vv[j];
            if (residual_degs[v] < best_deg) {
                best_deg = residual_degs[v];
                best_v_pos = j;
            }
        }

        int v = vv[best_v_pos];
        std::swap(vv[best_v_pos], vv[i]);

        for (int neighbour : g.adjlist[v])
            --residual_degs[neighbour];
    }

    vv.resize(new_sz);
}
auto simple_mwc(const SparseGraph & g, const Params params, VtxList & incumbent, long & search_node_count) -> void
{
    VtxList C(g.n);
    vector<int> vv(g.n);
    for (unsigned i=0; i<g.n; i++)
        vv[i] = i;
    if (params.unweighted_sort)
        reduce_and_reverse_unweighted_cp_order(vv, g);
    else
        reduce_and_reverse_cp_order(vv, g, 0, 0);
    if (vv.size()) {
        Graph ordered_graph = g.induced_subgraph<Graph>(vv);
        SimpleMWC(ordered_graph, params, incumbent).run(C, search_node_count);
        // Use vertex indices from original graph
        for (unsigned i=0; i<incumbent.vv.size(); i++)
            incumbent.vv[i] = vv[incumbent.vv[i]];
    }
}


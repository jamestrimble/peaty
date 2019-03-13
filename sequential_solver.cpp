#include "bitset.h"
#include "colourer.h"
#include "sequential_solver.h"
#include "graph.h"
#include "util.h"
#include "root_node_processing.h"
#include "params.h"

#include <algorithm>
#include <memory>
#include <vector>

auto update_incumbent_if_necessary(VtxList & C, VtxList & incumbent, const Params & params)
{
    if (C.total_wt > incumbent.total_wt) {
        incumbent = C;
//        if (!params.quiet)
//            printf("New incumbent of weight %ld\n", incumbent.total_wt);
    }
}

class MWC {
    Graph & g;
    const Params params;
    Colourer & colourer;
    VtxList & incumbent;
    vector<vector<unsigned long long>> branch_vv_bitsets;
    vector<vector<unsigned long long>> new_P_bitsets;

    void expand(VtxList& C, vector<unsigned long long> & P_bitset, long & search_node_count,
            std::atomic_int & ind_set_upper_bound)
    {
        ++search_node_count;

        if (bitset_empty(P_bitset, g.numwords)) {
            update_incumbent_if_necessary(C, incumbent, params);
            return;
        }

        int ind_set_upper_bound_int = ind_set_upper_bound.load();
        if (ind_set_upper_bound_int != -1 && int(incumbent.vv.size()) == ind_set_upper_bound_int) {
            return;
        }

        vector<unsigned long long> & branch_vv_bitset = branch_vv_bitsets[C.vv.size() - 1];
        if (branch_vv_bitset.empty())
            branch_vv_bitset.resize(g.numwords);
        else
            std::fill(branch_vv_bitset.begin(), branch_vv_bitset.end(), 0);

        long target = incumbent.total_wt - C.total_wt;
        if (colourer.colouring_bound(P_bitset, branch_vv_bitset, target)) {
            vector<unsigned long long> & new_P_bitset = new_P_bitsets[C.vv.size() - 1];
            if (new_P_bitset.empty())
                new_P_bitset.resize(g.numwords);

            bitset_intersect_with_complement(P_bitset, branch_vv_bitset, g.numwords);

            int v;
            while ((v=first_set_bit(branch_vv_bitset, g.numwords))!=-1) {
                unset_bit(branch_vv_bitset, v);
                bitset_intersection_with_complement(P_bitset, g.bit_complement_nd[v], new_P_bitset, g.numwords);
                C.push_vtx(v, g);
                expand(C, new_P_bitset, search_node_count, ind_set_upper_bound);
                set_bit(P_bitset, v);
                C.pop_vtx(g);
            }
        }
    }

public:
    MWC(Graph & g, const Params params, VtxList & incumbent, Colourer & colourer)
            : g(g), params(params), colourer(colourer), incumbent(incumbent), branch_vv_bitsets(g.n), new_P_bitsets(g.n)
    {
    }

    auto run(VtxList & C, long & search_node_count, std::atomic_int & ind_set_upper_bound) -> void
    {
        vector<unsigned long long> P_bitset(g.numwords, 0);
        set_first_n_bits(P_bitset, g.n);
        expand(C, P_bitset, search_node_count, ind_set_upper_bound);
    }
};

auto sequential_mwc(const SparseGraph & g, const Params params, VtxList & incumbent, long & search_node_count,
        std::atomic_int & ind_set_upper_bound) -> void
{
    VtxList C(g.n);

    auto vv0 = initialise(g);
//    printf("Initial incumbent weight %ld\n", incumbent.total_wt);
    remove_vertices_with_closed_nd_wt_leq_incumbent(g, vv0, 0, incumbent.total_wt, 1);
    SparseGraph ordered_graph = g.induced_subgraph<SparseGraph>(vv0);
    ordered_graph.sort_adj_lists();

    ++search_node_count;

    for (unsigned v=0; v<ordered_graph.n; v++) {
        C.push_vtx(v, ordered_graph.weight[v]);
        auto & v_adjlist = ordered_graph.adjlist[v];
        // vv is the list of vertices non-adjacent to v and less than v
        vector<int> vv;
        vector<bool> vv_bools(ordered_graph.n, true);
        for (int w : v_adjlist) {
            vv_bools[w] = false;
        }
        for (unsigned i=0; i<v; i++) {
            if (vv_bools[i]) {
                vv.push_back(i);
            }
        }
        long prev_incumbent_weight = incumbent.total_wt;
        update_incumbent_if_necessary(C, incumbent, params);
        // bear in mind that the adjacency lists of unordered_subgraph are not sorted!
        SparseGraph unordered_subgraph =
                ordered_graph.induced_subgraph<SparseGraph>(vv);
//        vector<int> vv1;
//        vv1.reserve(unordered_subgraph.n);
//        for (unsigned i=0; i<unordered_subgraph.n; i++)
//            vv1.push_back(i);
        auto vv1 = initialise(unordered_subgraph);
        if (vv1.size()) {
            Graph ordered_subgraph = unordered_subgraph.complement_of_induced_subgraph(vv1);
            std::shared_ptr<Colourer> colourer = Colourer::create_colourer(ordered_subgraph, params);
            MWC(ordered_subgraph, params, incumbent, *colourer).run(C, search_node_count, ind_set_upper_bound);
        }
        if (incumbent.total_wt > prev_incumbent_weight) {
            // Use vertex indices from original graph
            incumbent.vv[0] = vv0[incumbent.vv[0]];
            for (unsigned i=1; i<incumbent.vv.size(); i++)
                incumbent.vv[i] = vv0[vv[vv1[incumbent.vv[i]]]];
        }
        C.pop_vtx(ordered_graph.weight[v]);
    }
//    vector<int> vv1;
//    vv1.reserve(ordered_graph.n);
//    for (unsigned i=0; i<ordered_graph.n; i++)
//        vv1.push_back(i);
//    Graph ordered_subgraph = ordered_graph.complement_of_induced_subgraph(vv1);
//    std::shared_ptr<Colourer> colourer = Colourer::create_colourer(ordered_subgraph, params);
//
//    MWC(ordered_subgraph, params, incumbent, *colourer).run(C, search_node_count, ind_set_upper_bound);
//    for (unsigned i=0; i<incumbent.vv.size(); i++)
//        incumbent.vv[i] = vv0[incumbent.vv[i]];
}

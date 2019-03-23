#include "bitset.h"
#include "colourer.h"
#include "sequential_solver.h"
#include "graph.h"
#include "util.h"
#include "root_node_processing.h"
#include "params.h"

#include <iostream>  // TODO: don't include this

#include <algorithm>
#include <memory>
#include <random>
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

        vector<unsigned long long> & branch_vv_bitset = branch_vv_bitsets[C.vv.size()];
        if (branch_vv_bitset.empty())
            branch_vv_bitset.resize(g.numwords);
        else
            std::fill(branch_vv_bitset.begin(), branch_vv_bitset.end(), 0);

        long target = incumbent.total_wt - C.total_wt;
        if (colourer.colouring_bound(P_bitset, branch_vv_bitset, target)) {
            vector<unsigned long long> & new_P_bitset = new_P_bitsets[C.vv.size()];
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

class FastSet
{
    vector<bool> in_set;
    vector<int> position_in_elements_list;
public:
    vector<int> elements;
    FastSet(int capacity) : in_set(capacity), position_in_elements_list(capacity) {}

    void add(int x)
    {
        if (in_set[x])
            return;
        in_set[x] = true;
        position_in_elements_list[x] = elements.size();
        elements.push_back(x);
    }

    void remove(int x)
    {
        if (!in_set[x])
            return;
        int pos = position_in_elements_list[x];
        std::swap(elements.back(), elements[pos]);
        position_in_elements_list[elements[pos]] = pos;
        elements.pop_back();
        in_set[x] = false;
    }

    bool contains(int x)
    {
        return in_set[x];
    }

    unsigned size()
    {
        return elements.size();
    }
};

class LocalSearcher
{
    const SparseGraph & g;
    vector<int> num_conflicts;
    FastSet set_of_vv_with_no_conflicts;
    FastSet set_of_vv_with_one_conflict;
    vector<bool> ind_set;
    int ind_set_size;
    int tabu_duration;
    int time;
    vector<int> last_time_changed;
    vector<int> incumbent;
    std::mt19937 mt19937;

public:
    LocalSearcher(const SparseGraph & g) : g(g), num_conflicts(g.n),
            set_of_vv_with_no_conflicts(g.n), set_of_vv_with_one_conflict(g.n), ind_set(g.n),
            tabu_duration(10), time(tabu_duration + 1), last_time_changed(g.n)
    {
        for (unsigned i=0; i<g.n; i++) {
            set_of_vv_with_no_conflicts.add(i);
            set_of_vv_with_one_conflict.remove(i);
        }
    }

    void reset()
    {
        std::fill(num_conflicts.begin(), num_conflicts.end(), 0);
        std::fill(ind_set.begin(), ind_set.end(), false);
        std::fill(last_time_changed.begin(), last_time_changed.end(), 0);
        ind_set_size = 0;
        for (unsigned i=0; i<g.n; i++) {
            set_of_vv_with_no_conflicts.add(i);
            set_of_vv_with_one_conflict.remove(i);
        }
    }

    void add_to_ind_set(int v)
    {
        ind_set[v] = true;
        ++ind_set_size;
        for (int w : g.adjlist[v]) {
            if (num_conflicts[w] == 0) {
                set_of_vv_with_no_conflicts.remove(w);
                set_of_vv_with_one_conflict.add(w);
            } else if (num_conflicts[w] == 1) {
                set_of_vv_with_one_conflict.remove(w);
            }
            ++num_conflicts[w];
        }
    }

    void remove_from_ind_set(int v)
    {
        last_time_changed[v] = time;
        ind_set[v] = false;
        --ind_set_size;
        for (int w : g.adjlist[v]) {
            --num_conflicts[w];
            if (num_conflicts[w] == 0) {
                set_of_vv_with_no_conflicts.add(w);
                set_of_vv_with_one_conflict.remove(w);
            } else if (num_conflicts[w] == 1) {
                set_of_vv_with_one_conflict.add(w);
            }
        }
    }

    bool permitted_by_tabu_rule(int v)
    {
        return time > last_time_changed[v] + tabu_duration;
    }

    void greedily_add_to_is()
    {
        if (int(set_of_vv_with_no_conflicts.size()) == ind_set_size)
            return;

        vector<int> vertices_without_conflict;
        for (int v : set_of_vv_with_no_conflicts.elements)
            if (!ind_set[v])
                vertices_without_conflict.push_back(v);

        std::shuffle(vertices_without_conflict.begin(), vertices_without_conflict.end(), mt19937);

        for (int v : vertices_without_conflict)
            if (num_conflicts[v] == 0 && (ind_set_size>=int(incumbent.size()) || permitted_by_tabu_rule(v)))
                add_to_ind_set(v);

        if (ind_set_size > int(incumbent.size())) {
            incumbent.clear();
            for (unsigned i=0; i<g.n; i++) {
                if (ind_set[i]) {
                    incumbent.push_back(i);
                }
            }
            std::cout << "incumbent " << incumbent.size() << std::endl;
        }
    }

    void do_swap_or_deletion()
    {
        vector<int> vertices_with_one_conflict;

        // occasionally, do a drop even if a swap is possible
        std::uniform_int_distribution<> distrib0(0, 20);
        if (distrib0(mt19937)) {
            for (int v : set_of_vv_with_one_conflict.elements)
                if (permitted_by_tabu_rule(v))
                    vertices_with_one_conflict.push_back(v);
        }

        if (vertices_with_one_conflict.empty()) {
            // do a drop
            vector<int> vertices_in_is;
            for (unsigned i=0; i<g.n; i++)
                if (ind_set[i])
                    vertices_in_is.push_back(i);
            if (!vertices_in_is.empty()) {
                std::uniform_int_distribution<> distrib(0, vertices_in_is.size() - 1);
                int v = vertices_in_is[distrib(mt19937)];
                remove_from_ind_set(v);
            }
        } else {
            // do a swap
            std::uniform_int_distribution<> distrib(0, vertices_with_one_conflict.size() - 1);
            int v = vertices_with_one_conflict[distrib(mt19937)];
            for (int w : g.adjlist[v]) {
                if (ind_set[w]) {
                    remove_from_ind_set(w);
                    add_to_ind_set(v);
                    break;
                }
            }
        }
    }

    int find_lower_bound()
    {
        int local_time_limit = 5000;
        int iter_num = 0;
        while (time < 10000000) {
            int local_time = 0;
            int local_best = 0;
            while (local_time < local_time_limit) {
                greedily_add_to_is();
                do_swap_or_deletion();
                do_swap_or_deletion();
                if (ind_set_size > local_best) {
                    local_best = ind_set_size;
                    local_time = 0;
                }
                ++local_time;
                ++time;
                if (0 == (time % 100000)) {
                    std::cout << "c time " << time << std::endl;
                }
            }
            local_time_limit = local_time_limit + local_time_limit / 1000;
            reset();
            ++iter_num;
        }

        return incumbent.size();
    }
};

auto sequential_mwc(const SparseGraph & g, const Params params, VtxList & incumbent, long & search_node_count,
        std::atomic_int & ind_set_upper_bound) -> void
{
    VtxList C(g.n);

    LocalSearcher ls(g);
    int lower_bound = ls.find_lower_bound();
    std::cout << "Lower bound " << lower_bound << std::endl;

    auto vv0 = initialise(g);
//    printf("Initial incumbent weight %ld\n", incumbent.total_wt);
    remove_vertices_with_closed_nd_wt_leq_incumbent(g, vv0, 0, incumbent.total_wt, 1);
    SparseGraph ordered_graph = g.induced_subgraph<SparseGraph>(vv0);
    ordered_graph.sort_adj_lists();

    vector<int> vv1;
    vv1.reserve(ordered_graph.n);
    for (unsigned i=0; i<ordered_graph.n; i++)
        vv1.push_back(i);
    Graph ordered_subgraph = ordered_graph.complement_of_induced_subgraph(vv1);
    std::shared_ptr<Colourer> colourer = Colourer::create_colourer(ordered_subgraph, params);

    vector<unsigned long long> branch_vv_bitset(ordered_subgraph.numwords); 
    vector<unsigned long long> P_bitset(ordered_subgraph.numwords); 
    set_first_n_bits(P_bitset, ordered_subgraph.n);
    int target = 1;
    while (colourer->colouring_bound(P_bitset, branch_vv_bitset, target))
        ++target;

    while (true) {
        std::cout << "target " << target << std::endl;
        incumbent.total_wt = target - 1;
        MWC(ordered_subgraph, params, incumbent, *colourer).run(C, search_node_count, ind_set_upper_bound);
        if (incumbent.total_wt == target)
            break;
        target--;
    }
    for (unsigned i=0; i<incumbent.vv.size(); i++)
        incumbent.vv[i] = vv0[incumbent.vv[i]];
}

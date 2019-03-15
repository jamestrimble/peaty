#ifndef COLOURER_H
#define COLOURER_H

#include "bitset.h"
#include "colourer.h"
#include "params.h"
#include "fast_int_queue.h"
#include "graph.h"
#include "int_stack_without_dups.h"
#include "util.h"

#include <algorithm>
#include <memory>
#include <vector>

using std::vector;
struct Clause {
    vector<int> vv;
    long weight;
    long remaining_wt;
    unsigned long long sorting_score;
};

struct ListOfClauses {
    vector<Clause> clause;
    int size;
    int capacity;
    ListOfClauses(int capacity)
            : clause(capacity), size(0), capacity(capacity) { }
    auto clear() -> void { size = 0; }
};

// Which clauses does each vertex belong to?
using ClauseMembership = vector<vector<int>>;

class UnitPropagator {
    const Graph & g;
    const Params params;

    FastIntQueue Q;
    IntStackWithoutDups I;
    IntStackWithoutDups iset;
    ClauseMembership cm;

    vector<int> vv_count;
    vector<int> remaining_vv_count;

    vector<unsigned char> vertex_has_been_propagated;  // really a vector of bools

    // in unit_propagate_once, every vertex has a clause index as its reason,
    // or -1 if the vertex does not have a reason
    vector<int> reason;

    int get_unique_remaining_vtx(const Clause & c);

    void create_inconsistent_set(int c_idx, ListOfClauses & cc);

    // Return whether an inconsistent set has been found
    auto propagate_vertex(ListOfClauses & cc, int v, int u_idx,
            const vector<unsigned long long> & P_bitset) -> bool;

    auto unit_propagate_once(ListOfClauses & cc, int first_clause_index, int first_v,
            const vector<unsigned long long> & P_bitset) -> void;

    void remove_from_clause_membership(int v, int clause_idx);

    long process_inconsistent_set(IntStackWithoutDups & iset, ListOfClauses & cc);

    unsigned get_max_clause_size(const ListOfClauses & cc);

public:
    UnitPropagator(Graph & g, const Params params)
            : g(g), params(params), Q(g.n), I(g.n), iset(g.n), cm(g.n), vv_count(g.n), remaining_vv_count(g.n),
              vertex_has_been_propagated(g.n), reason(g.n)
    {
    }

    long unit_propagate(ListOfClauses & cc, long target_reduction, const vector<unsigned long long> & P_bitset);

    int unit_propagate_m1(ListOfClauses & cc, long target_reduction, long target,
            const vector<unsigned long long> & P_bitset);
};

/*******************************************************************************
*******************************************************************************/

class Colourer {
public:
    virtual auto colouring_bound(vector<unsigned long long> & P_bitset,
            vector<unsigned long long> & branch_vv_bitset, long target) -> bool = 0;

    static std::shared_ptr<Colourer> create_colourer(Graph & g, const Params & params);
};

class ClassEnlargingUnitPropColourer : public Colourer {
    Graph & g;
    const Params params;
    UnitPropagator unit_propagator;

    ListOfClauses cc;

    vector<int> vv;
    vector<unsigned long long> to_colour;
    vector<vector<unsigned long long>> candidates;
    vector<long> residual_wt;

public:
    ClassEnlargingUnitPropColourer(Graph & g, const Params params)
            : g(g), params(params), unit_propagator(g, params), cc(g.n),
              to_colour(g.numwords),
              candidates(2, vector<unsigned long long>(g.numwords)),
              residual_wt(g.n)
    {
    }

    auto try_to_enlarge_clause(Clause & clause, int numwords,
            vector<unsigned long long> & candidates, vector<unsigned long long> & to_colour) -> void
    {
        vv.clear();

        bitset_foreach(candidates, [this](int v){ vv.push_back(v); }, numwords);

        int sz = vv.size();
        for (int sum=0; sum<=sz*2-3; sum++) {
            int i_start = sum - sz + 1;
            if (i_start < 0)
                i_start = 0;
            for (int i=i_start, j=sum-i_start; i<j; i++, j--) {
                int w = vv[i];
                int u = vv[j];
                if (test_bit(g.bit_complement_nd[w], u)) {
                    clause.vv.pop_back();
                    clause.vv.push_back(w);
                    clause.vv.push_back(u);
                    return;
                }
            }
        }
    }
    bool colouring_bound(vector<unsigned long long> & P_bitset,
            vector<unsigned long long> & branch_vv_bitset, long target)
    {
        int numwords = calc_numwords(P_bitset, g.numwords);

        copy_bitset(P_bitset, to_colour, numwords);
        residual_wt = g.weight;
        cc.clear();

        long bound = 0;
        int v;
        int w = 0;
        while ((v=first_set_bit(to_colour, numwords))!=-1) {
            Clause & clause = cc.clause[cc.size];
            clause.vv.clear();
            clause.vv.push_back(v);
            bitset_intersection(to_colour, g.bit_complement_nd[v], candidates[0], numwords);
            int i = 0;
            while ((v=first_set_bit(candidates[i], numwords))!=-1) {
                clause.vv.push_back(v);
                bitset_intersection(candidates[i], g.bit_complement_nd[v], candidates[!i], numwords);
                i = !i;
                w = v;
            }
            if (clause.vv.size() > 1) {
                unset_bit(candidates[!i], w);
                try_to_enlarge_clause(clause, numwords, candidates[!i], to_colour);
            }
            long class_min_wt = residual_wt[clause.vv[0]];
            for (unsigned i=1; i<clause.vv.size(); i++) {
                int w = clause.vv[i];
                if (residual_wt[w] < class_min_wt)
                    class_min_wt = residual_wt[w];
            }

            for (int w : clause.vv) {
                residual_wt[w] -= class_min_wt;
                unset_bit_if(to_colour, w, residual_wt[w] <= 0);
            }
            bound += class_min_wt;
            clause.weight = class_min_wt;
            cc.size++;
        }

        for (unsigned long long i=0; i<cc.size; i++) {
            auto & c = cc.clause[i];
            c.sorting_score = ((unsigned long long)c.vv.size() << 32) - i;
        }
        std::sort(cc.clause.begin(), std::next(cc.clause.begin(), cc.size),
                [](auto & a, auto & b){return a.sorting_score > b.sorting_score;});

        long improvement = unit_propagator.unit_propagate(cc, bound-target, P_bitset);

        bool proved_we_can_prune = bound-improvement <= target;

        if (!proved_we_can_prune) {
            bound = 0;
            for (int i=0; i<cc.size; i++) {
                Clause & clause = cc.clause[i];
                assert (clause.weight >= 0);
                bound += clause.weight;
                if (bound > target)
                    for (int w : clause.vv)
                        set_bit(branch_vv_bitset, w);
            }
        }
        return !proved_we_can_prune;
    }
};

class UnitPropColourer : public Colourer {
    Graph & g;
    const Params params;
    UnitPropagator unit_propagator;

    ListOfClauses cc;

    vector<unsigned long long> to_colour;
    vector<unsigned long long> candidates;
    vector<long> residual_wt;

public:
    UnitPropColourer(Graph & g, const Params params)
            : g(g), params(params), unit_propagator(g, params), cc(g.n),
              to_colour(g.numwords),
              candidates(g.numwords),
              residual_wt(g.n)
    {
    }

    bool colouring_bound(vector<unsigned long long> & P_bitset,
            vector<unsigned long long> & branch_vv_bitset, long target)
    {
        int numwords = calc_numwords(P_bitset, g.numwords);

        copy_bitset(P_bitset, to_colour, numwords);
        residual_wt = g.weight;
        cc.clear();

        long bound = 0;
        int v;
        while ((v=first_set_bit(to_colour, numwords))!=-1) {
            Clause & clause = cc.clause[cc.size];
            clause.vv.clear();
            clause.vv.push_back(v);
            long class_min_wt = residual_wt[v];
            bitset_intersection(to_colour, g.bit_complement_nd[v], candidates, numwords);
            while ((v=first_set_bit(candidates, numwords))!=-1) {
                if (residual_wt[v] < class_min_wt)
                    class_min_wt = residual_wt[v];
                clause.vv.push_back(v);
                bitset_intersect_with(candidates, g.bit_complement_nd[v], numwords);
            }

            for (int w : clause.vv) {
                residual_wt[w] -= class_min_wt;
                unset_bit_if(to_colour, w, residual_wt[w] == 0);
            }
            bound += class_min_wt;
            clause.weight = class_min_wt;
            cc.size++;
        }

        long improvement = unit_propagator.unit_propagate(cc, bound-target, P_bitset);

        bool proved_we_can_prune = bound-improvement <= target;

        if (!proved_we_can_prune) {
            bound = 0;
            for (int i=0; i<cc.size; i++) {
                Clause & clause = cc.clause[i];
                assert (clause.weight >= 0);
                bound += clause.weight;
                if (bound > target)
                    for (int w : clause.vv)
                        set_bit(branch_vv_bitset, w);
            }
        }
        return !proved_we_can_prune;
    }
};

class EvenSimplerColourer : public Colourer {
    Graph & g;
    const Params params;

    vector<unsigned long long> candidates;
    vector<long> residual_wt;
    vector<int> col_class;

public:
    EvenSimplerColourer(Graph & g, const Params params)
            : g(g), params(params),
              candidates(g.numwords),
              residual_wt(g.n)
    {
    }

    auto colouring_bound(vector<unsigned long long> & P_bitset,
            vector<unsigned long long> & branch_vv_bitset, long target) -> bool
    {
        int numwords = calc_numwords(P_bitset, g.numwords);

        copy_bitset(P_bitset, branch_vv_bitset, numwords);
        residual_wt = g.weight;

        long bound = 0;
        int v;
        while ((v=first_set_bit(branch_vv_bitset, numwords))!=-1) {
            long class_min_wt = residual_wt[v];
            col_class.clear();
            col_class.push_back(v);
            bitset_intersection(branch_vv_bitset, g.bit_complement_nd[v], candidates, numwords);
            while ((v=first_set_bit(candidates, numwords))!=-1) {
                if (residual_wt[v] < class_min_wt)
                    class_min_wt = residual_wt[v];
                col_class.push_back(v);
                bitset_intersect_with(candidates, g.bit_complement_nd[v], numwords);
            }
            bound += class_min_wt;
            if (bound > target) {
                return true;
            }
            for (int w : col_class) {
                residual_wt[w] -= class_min_wt;
                unset_bit_if(branch_vv_bitset, w, residual_wt[w] == 0);
            }
        }
        return false;
    }
};

#endif

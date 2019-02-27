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

#ifdef VERY_VERBOSE
        printf("VERY_VERBOSE {\"clauses\": [");
        const char *sep = "";
        long total_wt = 0;
        for (int i=0; i<cc.size; i++) {
            printf("%s", sep);
            sep = ", ";
            Clause & c = cc.clause[i];
            printf("{\"weight\": %ld, ", c.weight);
            printf("\"vertices\": [");
            total_wt += c.weight;
            const char *sep2 = "";
            for (int v : c.vv) {
                printf("%s%d", sep2, v);
                sep2 = ", ";
            }
            printf("]");
            printf("}");
        }
        printf("], \"total_wt\": %ld}\n", total_wt);
#endif

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

#ifdef VERY_VERBOSE
        printf("VERY_VERBOSE {\"clauses\": [");
        const char *sep = "";
        long total_wt = 0;
        for (int i=0; i<cc.size; i++) {
            printf("%s", sep);
            sep = ", ";
            Clause & c = cc.clause[i];
            printf("{\"weight\": %ld, ", c.weight);
            printf("\"vertices\": [");
            total_wt += c.weight;
            const char *sep2 = "";
            for (int v : c.vv) {
                printf("%s%d", sep2, v);
                sep2 = ", ";
            }
            printf("]");
            printf("}");
        }
        printf("], \"total_wt\": %ld}\n", total_wt);
#endif

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

class M1UnitPropColourer : public Colourer {
    Graph & g;
    const Params params;
    UnitPropagator unit_propagator;

    ListOfClauses cc;

    vector<unsigned long long> to_colour;
    vector<unsigned long long> candidates;
    vector<long> residual_wt;

public:
    M1UnitPropColourer(Graph & g, const Params params)
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

        int start = unit_propagator.unit_propagate_m1(cc, bound-target, target, P_bitset);
        // start is the first clause C such that we haven't proved that all clauses up
        // to and including C must have an upper bound <= target

        for (int i=start; i<cc.size; i++) {
            Clause & clause = cc.clause[i];
            assert (clause.weight >= 0);
            for (int w : clause.vv)
                set_bit(branch_vv_bitset, w);
        }

        return start != cc.size;
    }
};

class MinMaxM1UnitPropColourer : public Colourer {
    Graph & g;
    const Params params;
    UnitPropagator unit_propagator;

    ListOfClauses cc;

    vector<unsigned long long> to_colour;
    vector<unsigned long long> to_colour2;
    vector<unsigned long long> candidates;
    vector<long> residual_wt;

public:
    MinMaxM1UnitPropColourer(Graph & g, const Params params)
            : g(g), params(params), unit_propagator(g, params), cc(g.n),
              to_colour(g.numwords),
              to_colour2(g.numwords),
              candidates(g.numwords),
              residual_wt(g.n)
    {
    }

    bool colouring_bound(vector<unsigned long long> & P_bitset,
            vector<unsigned long long> & branch_vv_bitset, long target)
    {
        int numwords = calc_numwords(P_bitset, g.numwords);

        copy_bitset(P_bitset, to_colour, numwords);
        clear_bitset(to_colour2, numwords);
        residual_wt = g.weight;
        cc.clear();

        long bound = 0;
        int v;
        while ((v=first_set_bit(to_colour, numwords))!=-1) {
            Clause & clause = cc.clause[cc.size];
            clause.vv.clear();
            clause.vv.push_back(v);
            long class_min_wt = residual_wt[v];
            long class_max_wt = residual_wt[v];
            bitset_intersection(to_colour, g.bit_complement_nd[v], candidates, numwords);
            while ((v=first_set_bit(candidates, numwords))!=-1) {
                if (residual_wt[v] < class_min_wt)
                    class_min_wt = residual_wt[v];
                if (residual_wt[v] > class_max_wt)
                    class_max_wt = residual_wt[v];
                clause.vv.push_back(v);
                bitset_intersect_with(candidates, g.bit_complement_nd[v], numwords);
            }

            if (bound + class_max_wt <= target) {
                bound += class_min_wt;
                clause.weight = class_min_wt;
                cc.size++;
                for (int w : clause.vv) {
                    residual_wt[w] -= class_min_wt;
                    unset_bit_if(to_colour, w, residual_wt[w] == 0);
                }
            } else {
                bitset_foreach(to_colour, [this, bound, target, &branch_vv_bitset](int v) {
                    if (bound + residual_wt[v] > target) {
                        set_bit(to_colour2, v);
                        unset_bit(to_colour, v);
                    }
                }, numwords);
            }
        }

        while ((v=first_set_bit(to_colour2, numwords))!=-1) {
            Clause & clause = cc.clause[cc.size];
            clause.vv.clear();
            clause.vv.push_back(v);
            long class_min_wt = residual_wt[v];
            bitset_intersection(to_colour2, g.bit_complement_nd[v], candidates, numwords);
            while ((v=first_set_bit(candidates, numwords))!=-1) {
                if (residual_wt[v] < class_min_wt)
                    class_min_wt = residual_wt[v];
                clause.vv.push_back(v);
                bitset_intersect_with(candidates, g.bit_complement_nd[v], numwords);
            }

            for (int w : clause.vv) {
                residual_wt[w] -= class_min_wt;
                unset_bit_if(to_colour2, w, residual_wt[w] == 0);
            }
            bound += class_min_wt;
            clause.weight = class_min_wt;
            cc.size++;
        }

        int start = unit_propagator.unit_propagate_m1(cc, bound-target, target, P_bitset);
        // start is the first clause C such that we haven't proved that all clauses up
        // to and including C must have an upper bound <= target

        for (int i=start; i<cc.size; i++) {
            Clause & clause = cc.clause[i];
            assert (clause.weight >= 0);
            for (int w : clause.vv)
                set_bit(branch_vv_bitset, w);
        }

        return start != cc.size;
    }
};

class SimpleColourer : public Colourer {
    Graph & g;
    const Params params;

    vector<unsigned long long> to_colour;
    vector<unsigned long long> candidates;
    vector<long> residual_wt;
    vector<int> col_class;

public:
    SimpleColourer(Graph & g, const Params params)
            : g(g), params(params),
              to_colour(g.numwords),
              candidates(g.numwords),
              residual_wt(g.n)
    {
    }

    auto colouring_bound(vector<unsigned long long> & P_bitset,
            vector<unsigned long long> & branch_vv_bitset, long target) -> bool
    {
        int numwords = calc_numwords(P_bitset, g.numwords);

        copy_bitset(P_bitset, to_colour, numwords);
        residual_wt = g.weight;

        long bound = 0;
        int v;
        bool target_exceeded = false;
        while ((v=first_set_bit(to_colour, numwords))!=-1) {
            long class_min_wt = residual_wt[v];
            long class_max_wt = residual_wt[v];
            col_class.clear();
            col_class.push_back(v);
            bitset_intersection(to_colour, g.bit_complement_nd[v], candidates, numwords);
            while ((v=first_set_bit(candidates, numwords))!=-1) {
                if (residual_wt[v] < class_min_wt)
                    class_min_wt = residual_wt[v];
                if (residual_wt[v] > class_max_wt)
                    class_max_wt = residual_wt[v];
                col_class.push_back(v);
                bitset_intersect_with(candidates, g.bit_complement_nd[v], numwords);
            }
            if (bound + class_max_wt <= target) {
                bound += class_min_wt;
                for (int w : col_class) {
                    residual_wt[w] -= class_min_wt;
                    unset_bit_if(to_colour, w, residual_wt[w] == 0);
                }
            } else {
                target_exceeded = true;
                bitset_foreach(to_colour, [this, bound, target, &branch_vv_bitset](int v) {
                    if (bound + residual_wt[v] > target) {
                        set_bit(branch_vv_bitset, v);
                        unset_bit(to_colour, v);
                    }
                }, numwords);
            }
        }
        return target_exceeded;
    }
};

class SimpleColourerWordByWord : public Colourer {
    Graph & g;
    const Params params;

    vector<unsigned long long> to_colour;
    vector<unsigned long long> candidates;
    vector<long> residual_wt;
    vector<int> col_class;

    vector<vector<unsigned long long>> bit_complement_nd_transposed;

public:
    SimpleColourerWordByWord(Graph & g, const Params params)
            : g(g), params(params),
              to_colour(g.numwords),
              candidates(g.numwords),
              residual_wt(g.n),
              bit_complement_nd_transposed(g.numwords)
    {
        for (int i=0; i<g.numwords; i++) {
            bit_complement_nd_transposed[i].reserve(g.n);
            for (int j=0; j<g.n; j++) {
                bit_complement_nd_transposed[i].push_back(g.bit_complement_nd[j][i]);
            }
        }
    }

    auto colouring_bound(vector<unsigned long long> & P_bitset,
            vector<unsigned long long> & branch_vv_bitset, long target) -> bool
    {
        int numwords = calc_numwords(P_bitset, g.numwords);

        copy_bitset(P_bitset, to_colour, numwords);
        residual_wt = g.weight;

        long bound = 0;
        bool target_exceeded = false;
        for (int i=0; i<numwords; i++) {
            unsigned long long & word = to_colour[i];
            while (word) {
                int bit = __builtin_ctzll(word);
                int u = i * BITS_PER_WORD + bit;
                long class_min_wt = residual_wt[u];
                long class_max_wt = residual_wt[u];
                col_class.clear();
                col_class.push_back(u);
                for (int j=i; j<numwords; j++) {
                    unsigned long long word2 = to_colour[j] & bit_complement_nd_transposed[j][u];
                    for (unsigned k=1; k<col_class.size(); k++) {
                        int w = col_class[k];
                        word2 &= bit_complement_nd_transposed[j][w];
                    }
                    while (word2) {
                        int bit2 = __builtin_ctzll(word2);
                        int v = j * BITS_PER_WORD + bit2;
                        if (residual_wt[v] < class_min_wt)
                            class_min_wt = residual_wt[v];
                        if (residual_wt[v] > class_max_wt)
                            class_max_wt = residual_wt[v];
                        col_class.push_back(v);
                        word2 &= bit_complement_nd_transposed[j][v];
                    }
                }
                if (bound + class_max_wt <= target) {
                    bound += class_min_wt;
                    for (int w : col_class) {
                        residual_wt[w] -= class_min_wt;
                        unset_bit_if(to_colour, w, residual_wt[w] == 0);
                    }
                } else {
                    target_exceeded = true;
                    bitset_foreach(to_colour, [this, bound, target, &branch_vv_bitset](int v) {
                        if (bound + residual_wt[v] > target) {
                            set_bit(branch_vv_bitset, v);
                            unset_bit(to_colour, v);
                        }
                    }, numwords);
                }
            }
        }
        return target_exceeded;
    }
};

class SimpleColourerWithSuperSimplePreCheck : public Colourer {
    Graph & g;
    const Params params;

    vector<unsigned long long> to_colour;
    vector<unsigned long long> candidates;
    vector<long> residual_wt;
    vector<int> col_class;

public:
    SimpleColourerWithSuperSimplePreCheck(Graph & g, const Params params)
            : g(g), params(params),
              to_colour(g.numwords),
              candidates(g.numwords),
              residual_wt(g.n)
    {
    }

    auto colouring_bound(vector<unsigned long long> & P_bitset,
            vector<unsigned long long> & branch_vv_bitset, long target) -> bool
    {
        int numwords = calc_numwords(P_bitset, g.numwords);

        copy_bitset(P_bitset, to_colour, numwords);

        long bound = 0;
        int v;
        bool target_exceeded = false;
        while ((v=first_set_bit(to_colour, numwords))!=-1) {
            long class_max_wt = g.weight[v];
            bitset_intersection(to_colour, g.bit_complement_nd[v], candidates, numwords);
            unset_bit(to_colour, v);
            while ((v=first_set_bit(candidates, numwords))!=-1) {
                if (g.weight[v] > class_max_wt)
                    class_max_wt = g.weight[v];
                unset_bit(to_colour, v);
                bitset_intersect_with(candidates, g.bit_complement_nd[v], numwords);
            }
            bound += class_max_wt;
            if (bound > target) {
                target_exceeded = true;
                break;
            }
        }

        if (!target_exceeded)
            return false;

        copy_bitset(P_bitset, to_colour, numwords);
        residual_wt = g.weight;

        bound = 0;
        target_exceeded = false;
        while ((v=first_set_bit(to_colour, numwords))!=-1) {
            long class_min_wt = residual_wt[v];
            long class_max_wt = residual_wt[v];
            col_class.clear();
            col_class.push_back(v);
            bitset_intersection(to_colour, g.bit_complement_nd[v], candidates, numwords);
            while ((v=first_set_bit(candidates, numwords))!=-1) {
                if (residual_wt[v] < class_min_wt)
                    class_min_wt = residual_wt[v];
                if (residual_wt[v] > class_max_wt)
                    class_max_wt = residual_wt[v];
                col_class.push_back(v);
                bitset_intersect_with(candidates, g.bit_complement_nd[v], numwords);
            }
            if (bound + class_max_wt <= target) {
                bound += class_min_wt;
                for (int w : col_class) {
                    residual_wt[w] -= class_min_wt;
                    unset_bit_if(to_colour, w, residual_wt[w] == 0);
                }
            } else {
                target_exceeded = true;
                bitset_foreach(to_colour, [this, bound, target, &branch_vv_bitset](int v) {
                    if (bound + residual_wt[v] > target) {
                        set_bit(branch_vv_bitset, v);
                        unset_bit(to_colour, v);
                    }
                }, numwords);
            }
        }
        return target_exceeded;
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

// linked-list clause (a clause that can appear in a doubly-linked list)
struct LLClauseUsingBitsets {
    vector<unsigned long long> bitset;
    long weight;
    LLClauseUsingBitsets * next;
    LLClauseUsingBitsets * prev;
};
struct LLClauseUsingVtxLists {
    vector<int> vv;
    long weight;
    LLClauseUsingVtxLists * next;
    LLClauseUsingVtxLists * prev;
};

template<typename T>
class ListOfLLClauses {
public:
    vector<T> clause;
    int size;
    int capacity;
    ListOfLLClauses(int capacity)
            : clause(capacity), size(0), capacity(capacity) { }
    auto clear() -> void { size = 0; }
};

class VertexByVertexColourerUsingBitsets : public Colourer {
    Graph & g;
    const Params params;

    ListOfLLClauses<LLClauseUsingBitsets> cc;

    vector<unsigned long long> to_colour;
    vector<unsigned long long> candidates;
    vector<long> residual_wt;
    vector<LLClauseUsingBitsets *> clause_ptrs;

public:
    VertexByVertexColourerUsingBitsets(Graph & g, const Params params)
            : g(g), params(params), cc(g.n),
              to_colour(g.numwords),
              residual_wt(g.n)
    {
    }

    void remove_clause_if_intersection_with_to_colour_is_empty(LLClauseUsingBitsets & clause, int numwords)
    {
        if (have_empty_intersection(clause.bitset, to_colour, numwords)) {
            clause.prev->next = clause.next;
            clause.next->prev = clause.prev;
        }
    }

    bool try_to_use_vtx(int w, LLClauseUsingBitsets & list_head, LLClauseUsingBitsets & list_tail,
            long & bound, long target, int numwords)
    {
        clause_ptrs.clear();
        long remaining_wt_to_cover = g.weight[w];
        for (LLClauseUsingBitsets * clause_ptr=list_head.next; clause_ptr!=&list_tail; clause_ptr=clause_ptr->next) {
            if (test_bit(clause_ptr->bitset, w)) {
                auto & clause = *clause_ptr;
                if (clause.weight < remaining_wt_to_cover) {
                    remaining_wt_to_cover -= clause.weight;
                    clause_ptrs.push_back(clause_ptr);
                } else {
                    if (clause.weight > remaining_wt_to_cover) {
                        auto & new_clause = cc.clause[cc.size];
                        new_clause.bitset = clause.bitset;
                        new_clause.weight = clause.weight - remaining_wt_to_cover;
                        new_clause.next = clause.next;
                        new_clause.prev = clause_ptr;
                        clause.weight = remaining_wt_to_cover;
                        clause.next->prev = &new_clause;
                        clause.next = &cc.clause[cc.size];
                        ++cc.size;
                    }
                    for (auto clause_ptr : clause_ptrs) {
                        bitset_intersect_with(clause_ptr->bitset, g.bit_complement_nd[w], numwords);
                        remove_clause_if_intersection_with_to_colour_is_empty(*clause_ptr, numwords);
                    }
                    bitset_intersect_with(clause.bitset, g.bit_complement_nd[w], numwords);
                    remove_clause_if_intersection_with_to_colour_is_empty(clause, numwords);
                    return true;
                }
            }
        }
        if (bound + remaining_wt_to_cover <= target) {
            for (auto clause_ptr : clause_ptrs) {
                bitset_intersect_with(clause_ptr->bitset, g.bit_complement_nd[w], numwords);
                remove_clause_if_intersection_with_to_colour_is_empty(*clause_ptr, numwords);
            }
            if (have_non_empty_intersection(g.bit_complement_nd[w], to_colour, numwords)) {
                LLClauseUsingBitsets & new_clause = cc.clause[cc.size];
                new_clause.bitset = g.bit_complement_nd[w];
                new_clause.weight = remaining_wt_to_cover;
                new_clause.next = &list_tail;
                new_clause.prev = list_tail.prev;
                list_tail.prev->next = &new_clause;
                list_tail.prev = &new_clause;
                ++cc.size;
            }
            bound += remaining_wt_to_cover;
            return true;
        }
        return false;
    }

    bool colouring_bound(vector<unsigned long long> & P_bitset,
            vector<unsigned long long> & branch_vv_bitset, long target)
    {
        int numwords = calc_numwords(P_bitset, g.numwords);

        if (numwords == 0)
            return false;

        copy_bitset(P_bitset, to_colour, numwords);

        residual_wt = g.weight;
        cc.clear();

        long bound = 0;

        LLClauseUsingBitsets list_head;
        LLClauseUsingBitsets list_tail;
        list_head.next = &list_tail;
        list_tail.prev = &list_head;

        bool target_exceeded = false;
        for (int i=0; i<numwords; i++) {
            unsigned long long & word = to_colour[i];
            while (word) {
                int bit = __builtin_ctzll(word);
                word ^= (1ull << bit);
                int w = i*BITS_PER_WORD + bit;
                if (!try_to_use_vtx(w, list_head, list_tail, bound, target, numwords)) {
                    set_bit(branch_vv_bitset, w);
                    target_exceeded = true;
                }
            }
        }

        return target_exceeded;
    }
};

class VertexByVertexColourer : public Colourer {
    Graph & g;
    const Params params;

    ListOfLLClauses<LLClauseUsingVtxLists> cc;

    vector<unsigned long long> to_colour;
    vector<unsigned long long> candidates;
    vector<long> residual_wt;
    vector<LLClauseUsingVtxLists *> clause_ptrs;

public:
    VertexByVertexColourer(Graph & g, const Params params)
            : g(g), params(params), cc(g.n),
              to_colour(g.numwords),
              residual_wt(g.n)
    {
    }

    bool fits(int w, const LLClauseUsingVtxLists & clause)
    {
        for (int u : clause.vv) {
            if (!test_bit(g.bit_complement_nd[w], u)) {
                return false;
            }
        }
        return true;
    }

    bool try_to_use_vtx(int w, LLClauseUsingVtxLists & list_head, LLClauseUsingVtxLists & list_tail,
            long & bound, long target, int numwords)
    {
        clause_ptrs.clear();
        long remaining_wt_to_cover = g.weight[w];
        for (LLClauseUsingVtxLists * clause_ptr=list_head.next; clause_ptr!=&list_tail; clause_ptr=clause_ptr->next) {
            if (fits(w, *clause_ptr)) {
                auto & clause = *clause_ptr;
                if (clause.weight < remaining_wt_to_cover) {
                    remaining_wt_to_cover -= clause.weight;
                    clause_ptrs.push_back(clause_ptr);
                } else {
                    if (clause.weight > remaining_wt_to_cover) {
                        auto & new_clause = cc.clause[cc.size];
                        new_clause.vv = clause.vv;
                        new_clause.weight = clause.weight - remaining_wt_to_cover;
                        new_clause.next = clause.next;
                        new_clause.prev = clause_ptr;
                        clause.weight = remaining_wt_to_cover;
                        clause.next->prev = &new_clause;
                        clause.next = &cc.clause[cc.size];
                        ++cc.size;
                    }
                    for (auto clause_ptr : clause_ptrs) {
                        clause_ptr->vv.push_back(w);
                    }
                    clause.vv.push_back(w);
                    return true;
                }
            }
        }
        if (bound + remaining_wt_to_cover <= target) {
            for (auto clause_ptr : clause_ptrs) {
                clause_ptr->vv.push_back(w);
            }
            if (have_non_empty_intersection(g.bit_complement_nd[w], to_colour, numwords)) {
                LLClauseUsingVtxLists & new_clause = cc.clause[cc.size];
                new_clause.vv.clear();
                new_clause.vv.push_back(w);
                new_clause.weight = remaining_wt_to_cover;
                new_clause.next = &list_tail;
                new_clause.prev = list_tail.prev;
                list_tail.prev->next = &new_clause;
                list_tail.prev = &new_clause;
                ++cc.size;
            }
            bound += remaining_wt_to_cover;
            return true;
        }
        return false;
    }

    bool colouring_bound(vector<unsigned long long> & P_bitset,
            vector<unsigned long long> & branch_vv_bitset, long target)
    {
        int numwords = calc_numwords(P_bitset, g.numwords);

        if (numwords == 0)
            return false;

        copy_bitset(P_bitset, to_colour, numwords);

        residual_wt = g.weight;
        cc.clear();

        long bound = 0;

        LLClauseUsingVtxLists list_head;
        LLClauseUsingVtxLists list_tail;
        list_head.next = &list_tail;
        list_tail.prev = &list_head;

        bool target_exceeded = false;
        for (int i=0; i<numwords; i++) {
            unsigned long long & word = to_colour[i];
            while (word) {
                int bit = __builtin_ctzll(word);
                word ^= (1ull << bit);
                int w = i*BITS_PER_WORD + bit;
                if (!try_to_use_vtx(w, list_head, list_tail, bound, target, numwords)) {
                    set_bit(branch_vv_bitset, w);
                    target_exceeded = true;
                }
            }
        }

        return target_exceeded;
    }
};

class HebrardKatsirelosColourer : public Colourer {
    Graph & g;
    const Params params;

    vector<unsigned long long> candidates;
    vector<unsigned long long> to_colour_this_layer;
    vector<long> residual_wt;
    vector<int> col_class;

public:
    HebrardKatsirelosColourer(Graph & g, const Params params)
            : g(g), params(params),
              candidates(g.numwords),
              to_colour_this_layer(g.numwords),
              residual_wt(g.n)
    {
    }

    // Note: this could be optimised a bit more.
    auto colouring_bound(vector<unsigned long long> & P_bitset,
            vector<unsigned long long> & branch_vv_bitset, long target) -> bool
    {
        int numwords = calc_numwords(P_bitset, g.numwords);

        copy_bitset(P_bitset, branch_vv_bitset, numwords);
        residual_wt = g.weight;

        long bound = 0;
        int v;

        bool allow_unit_classes = false;

        while (!bitset_empty(branch_vv_bitset, numwords)) {
            copy_bitset(branch_vv_bitset, to_colour_this_layer, numwords);
            bool found_non_unit_class = false;
            while ((v=first_set_bit(to_colour_this_layer, numwords))!=-1) {
                unset_bit(to_colour_this_layer, v);
                long class_min_wt = residual_wt[v];
                col_class.clear();
                col_class.push_back(v);
                bitset_intersection(branch_vv_bitset, g.bit_complement_nd[v], candidates, numwords);
                while ((v=first_set_bit(candidates, numwords))!=-1) {
                    unset_bit(candidates, v);
                    if (residual_wt[v] < class_min_wt)
                        class_min_wt = residual_wt[v];
                    col_class.push_back(v);
                    bitset_intersect_with(candidates, g.bit_complement_nd[v], numwords);
                }
                if (col_class.size() > 1 || allow_unit_classes) {
                    bound += class_min_wt;
                    if (bound > target) {
                        return true;
                    }
                    for (int w : col_class) {
                        residual_wt[w] -= class_min_wt;
                        unset_bit_if(branch_vv_bitset, w, residual_wt[w] == 0);
                    }
                    if (col_class.size() > 1)
                        found_non_unit_class = true;
                }
            }
            if (!found_non_unit_class)
                allow_unit_classes = true;
        }
        return false;
    }
};

#endif

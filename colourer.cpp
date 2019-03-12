#include "colourer.h"

std::shared_ptr<Colourer> Colourer::create_colourer(Graph & g, const Params & params)
{
    switch(params.colouring_variant)
    {
    case 2:
        return std::make_shared<UnitPropColourer>(g, params);
    case 3:
        return std::make_shared<ClassEnlargingUnitPropColourer>(g, params);
    default:   // option 0
        return std::make_shared<EvenSimplerColourer>(g, params);
    }
}

int UnitPropagator::get_unique_remaining_vtx(const Clause & c) {
    unsigned i = 0;
    int v;
    while (reason[v = c.vv[i]] != -1)
        ++i;
    assert(i < c.vv.size());
    return v;
}

void UnitPropagator::create_inconsistent_set(int c_idx, ListOfClauses & cc)
{
    I.push(c_idx);
    unsigned j = 0;
    while (j != I.vals.size()) {
        for (int w : cc.clause[I.vals[j]].vv) {
            int r = reason[w];
            if (r != -1 && !I.on_stack[r]) {
                I.push(r);
            }
        }
        ++j;
    }
}

// Return whether an inconsistent set has been found
auto UnitPropagator::propagate_vertex(ListOfClauses & cc, int v, int u_idx,
        const vector<unsigned long long> & P_bitset) -> bool
{
    for (int i=0; i<g.numwords; i++) {
        // iterate over vertices w that are in both P and the complement-graph
        // neighbourhood of v
        unsigned long long word = g.bit_complement_nd[v][i] & P_bitset[i];
        while (word) {
            int bit = __builtin_ctzll(word);
            word ^= (1ull << bit);
            int w = i*BITS_PER_WORD + bit;
            if (reason[w] == -1) {
                reason[w] = u_idx;
                for (int c_idx : cm[w]) {
                    remaining_vv_count[c_idx]--;
                    if (remaining_vv_count[c_idx] == 1) {
                        Q.enqueue(c_idx);
                    } else if (remaining_vv_count[c_idx] == 0) {
                        create_inconsistent_set(c_idx, cc);
                        return true;
                    }
                }
            }
        }
    }
    return false;
}

auto UnitPropagator::unit_propagate_once(ListOfClauses & cc, int first_clause_index, int first_v,
        const vector<unsigned long long> & P_bitset) -> void
{
    I.clear();
    Q.clear();

    remaining_vv_count = vv_count;
    std::fill(vertex_has_been_propagated.begin(), vertex_has_been_propagated.end(), 0);
    std::fill(reason.begin(), reason.end(), -1);

    if (propagate_vertex(cc, first_v, first_clause_index, P_bitset))
        return;
    vertex_has_been_propagated[first_v] = true;

    while (!Q.empty()) {
        int u_idx = Q.dequeue();
        assert (remaining_vv_count[u_idx] == 1);
        int v = get_unique_remaining_vtx(cc.clause[u_idx]);
        if (!vertex_has_been_propagated[v]) {
            if (propagate_vertex(cc, v, u_idx, P_bitset))
                return;
            vertex_has_been_propagated[v] = true;
        }
    }
}

void UnitPropagator::remove_from_clause_membership(int v, int clause_idx)
{
    assert(std::find(cm[v].begin(), cm[v].end(), clause_idx) != cm[v].end());
    cm[v].erase(std::find(cm[v].begin(), cm[v].end(), clause_idx));
}

long UnitPropagator::process_inconsistent_set(IntStackWithoutDups & iset, ListOfClauses & cc)
{
    assert(!iset.vals.empty());

    // find max index and min remaining weight
    int max_idx = iset.vals[0];
    long min_wt = cc.clause[max_idx].remaining_wt;
    for (unsigned i=1; i<iset.vals.size(); i++) {
        int c_idx = iset.vals[i];
        long wt = cc.clause[c_idx].remaining_wt;
        if (wt < min_wt)
            min_wt = wt;
        if (c_idx > max_idx)
            max_idx = c_idx;
    }

    for (int c_idx : iset.vals) {
        Clause & c = cc.clause[c_idx];
        c.remaining_wt -= min_wt;
        if (c.remaining_wt == 0) {
            // Remove references to this clause from CM
            for (int v : c.vv)
                remove_from_clause_membership(v, c_idx);
        }
    }
    cc.clause[max_idx].weight -= min_wt;  // decrease weight of last clause in set
    return min_wt;
}

unsigned UnitPropagator::get_max_clause_size(const ListOfClauses & cc)
{
    unsigned max_size = 0;
    for (int i=0; i<cc.size; i++) {
        unsigned sz = cc.clause[i].vv.size();
        if (sz > max_size)
            max_size = sz;
    }
    return max_size;
}

long UnitPropagator::unit_propagate(ListOfClauses & cc, long target_reduction, const vector<unsigned long long> & P_bitset)
{
    if (target_reduction <= 0)
        return 0;

    for (int v=0; v<g.n; v++)
        cm[v].clear();

    for (int i=0; i<cc.size; i++) {
        Clause & clause = cc.clause[i];
        vv_count[i] = clause.vv.size();
        for (int v : clause.vv)
            cm[v].push_back(i);
    }

    for (int i=0; i<cc.size; i++)
        cc.clause[i].remaining_wt = cc.clause[i].weight;

    long improvement = 0;

#ifdef VERY_VERBOSE
    printf("VERY_VERBOSE {\"isets\": [");
    const char *sep = "";
#endif
    unsigned max_clause_size = params.max_sat_level == -1 ?
            get_max_clause_size(cc) : params.max_sat_level;

    for (unsigned clause_size = 1; clause_size <= max_clause_size; clause_size++) {
        for (int i=0; i<cc.size; i++) {
            Clause & clause = cc.clause[i];
            if (clause.vv.size() != clause_size)
                continue;

            for (;;) {
                if (clause.remaining_wt == 0)
                    break;

                iset.clear();
                for (int v : clause.vv) {
                    unit_propagate_once(cc, i, v, P_bitset);
                    if (I.vals.empty()) {
                        iset.clear();
                        break;
                    }
                    for (int clause_idx : I.vals)
                        iset.push(clause_idx);
                }

                if (iset.empty())
                    break;

#ifdef VERY_VERBOSE
                printf("%s[", sep);
                sep = ", ";
                const char *sep2 = "";
                for (int val : iset.vals) {        // TODO: rename val to something more meaningful
                    printf("%s%d", sep2, val);
                    sep2 = ", ";
                }
                printf("]");
#endif
                improvement += process_inconsistent_set(iset, cc);

                if (improvement >= target_reduction)
                    return improvement;
            }
        }
    }
#ifdef VERY_VERBOSE
    printf("]}\n");
#endif

    return improvement;
}

int UnitPropagator::unit_propagate_m1(ListOfClauses & cc, long target_reduction, long target,
        const vector<unsigned long long> & P_bitset)
{
    if (target_reduction <= 0)
        return cc.size;

    for (int v=0; v<g.n; v++)
        cm[v].clear();

    for (int i=0; i<cc.size; i++) {
        Clause & clause = cc.clause[i];
        vv_count[i] = clause.vv.size();
        for (int v : clause.vv)
            cm[v].push_back(i);
    }

    for (int i=0; i<cc.size; i++)
        cc.clause[i].remaining_wt = cc.clause[i].weight;

    long improvement = 0;

    long bound = 0;
    for (int i=0; i<cc.size; i++) {
        Clause & clause = cc.clause[i];
        if (clause.vv.size() == 1) {
            for (;;) {
                if (clause.remaining_wt == 0)
                    break;

                int v = clause.vv[0];
                unit_propagate_once(cc, i, v, P_bitset);
                if (I.empty())
                    break;

                improvement += process_inconsistent_set(I, cc);

                if (improvement >= target_reduction)
                    return cc.size;
            }
        }

        bound += clause.weight;
        if (bound > target)
            return i;
    }

    return cc.size;
}

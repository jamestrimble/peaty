#include "bitset.h"
#include "sequential_solver.h"
#include "graph.h"
#include "util.h"

#include <algorithm>
#include <thread>
#include <vector>

vector<int> initialise(const SparseGraph & g) {
//    printf("Starting init\n");
    vector<int> vv;

    if (!g.n)
        return vv;

    vector<int> residual_degs;
    residual_degs.reserve(g.n);
    for (auto & lst : g.adjlist)
        residual_degs.push_back(lst.size());

    vector<int> ll_next(g.n * 2);
    vector<int> ll_prev(g.n * 2);
    int *list_head = &ll_next[g.n];
    int *list_tail = &ll_prev[g.n];
    for (unsigned i=0; i<g.n; i++) {
        list_head[i] = g.n + i;
        list_tail[i] = g.n + i;
    }
    for (unsigned i=0; i<g.n; i++) {
        int deg = residual_degs[i];
        ll_prev[i] = g.n + deg;
        ll_next[i] = list_head[deg];
        ll_prev[ll_next[i]] = i;
        list_head[deg] = i;
    }
//    printf("Made lists\n");

    unsigned list_idx = g.n - 1;

    vector<unsigned char> in_vv(g.n, 0);
    for (;;) {
        while (list_head[list_idx] >= int(g.n))
            --list_idx;

        if (list_idx == 0) {
            for (unsigned v=0; v<g.n; v++) {
                if (!in_vv[v]) {
                    vv.push_back(v);
                }
            }
            std::reverse(vv.begin(), vv.end());
            return vv;
        }

        int v = list_head[list_idx];
        vv.push_back(v);
        in_vv[v] = 1;

        // remove from list
        ll_next[ll_prev[v]] = ll_next[v];
        ll_prev[ll_next[v]] = ll_prev[v];

        for (int neighbour : g.adjlist[v]) {
            if (!in_vv[neighbour]) {
                // remove from list
                ll_next[ll_prev[neighbour]] = ll_next[neighbour];
                ll_prev[ll_next[neighbour]] = ll_prev[neighbour];

                --residual_degs[neighbour];
                int r = residual_degs[neighbour];

                // insert into list
                ll_prev[neighbour] = g.n + r;
                ll_next[neighbour] = list_head[r];
                ll_prev[ll_next[neighbour]] = neighbour;
                list_head[r] = neighbour;
            }
        }
        if (list_idx < g.n - 1)
            ++list_idx;
    }
}

vector<long> calc_weighted_degs(const SparseGraph & g) {
    vector<long> weighted_degs(g.n);
    for (int v=0; v<int(g.n); v++)
        for (auto w : g.adjlist[v])
            weighted_degs[w] += g.weight[v];

    return weighted_degs;
}

auto remove_vertices_with_closed_nd_wt_leq_incumbent(const SparseGraph & g,
        vector<int> & vv, long current_wt, long incumbent_wt, int num_threads) -> void
{
    vector<long> weighted_degs = calc_weighted_degs(g);

    vv.erase(std::remove_if(vv.begin(), vv.end(),
            [&g, &weighted_degs, current_wt, incumbent_wt](int v) {
                return g.weight[v] + weighted_degs[v] + current_wt <= incumbent_wt;
            }), vv.end());
}

auto reduce_and_reverse_cp_order(vector<int> & vv, const SparseGraph & g,
        long current_wt, long incumbent_wt) -> void
{
    if (vv.empty())
        return;

    auto residual_weighted_degs = calc_weighted_degs(g);

    unsigned new_sz = vv.size();

    for (unsigned i=vv.size(); i--; ) {
        // find vertex with lowest residual_weighted_deg
        unsigned best_v_pos = i;
        long best_wt_deg = residual_weighted_degs[vv[i]];
        for (unsigned j=i; j--; ) {
            int v = vv[j];
            if (residual_weighted_degs[v] < best_wt_deg) {
                best_wt_deg = residual_weighted_degs[v];
                best_v_pos = j;
            }
        }

        int v = vv[best_v_pos];
        std::swap(vv[best_v_pos], vv[i]);

        if (new_sz==i+1 && g.weight[v] + best_wt_deg + current_wt <= incumbent_wt)
            --new_sz;

        for (int neighbour : g.adjlist[v])
            residual_weighted_degs[neighbour] -= g.weight[v];
    }

    vv.resize(new_sz);
}



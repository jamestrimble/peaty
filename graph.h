#ifndef GRAPH_H
#define GRAPH_H

#include <limits.h>
#include <stdbool.h>
#include <stdlib.h>

#include <vector>

using std::vector;

#define BYTES_PER_WORD sizeof(unsigned long long)
#define BITS_PER_WORD (CHAR_BIT * BYTES_PER_WORD)

#define SET_BIT(bitset, bit) do { bitset[bit/BITS_PER_WORD] |= (1ull << (bit%BITS_PER_WORD)); } while (0);
#define UNSET_BIT(bitset, bit) do { bitset[bit/BITS_PER_WORD] &= ~(1ull << (bit%BITS_PER_WORD)); } while (0);

struct Graph
{
    int n;
    int numwords;
    vector<long> weight;
    std::vector<std::vector<unsigned long long>> bit_complement_nd;

    Graph(int n);
    auto add_edge(int v, int w) -> void
    {
        UNSET_BIT(bit_complement_nd[v], w);
        UNSET_BIT(bit_complement_nd[w], v);
    }

    auto remove_edge(int v, int w) -> void
    {
        SET_BIT(bit_complement_nd[v], w);
        SET_BIT(bit_complement_nd[w], v);
    }

    auto resize(unsigned new_n) -> void;
};

struct VtxList {
    long total_wt;
    std::vector<int> vv;
    VtxList(int capacity) : total_wt(0)
    {
        vv.reserve(capacity);
    }
    auto clear() -> void
    {
        total_wt = 0;
        vv.clear();
    }
    auto push_vtx(int v, const Graph & g) -> void
    {
        vv.push_back(v);
        total_wt += g.weight[v];
    }
    auto pop_vtx(const Graph & g) -> void
    {
        total_wt -= g.weight[vv.back()];
        vv.pop_back();
    }
    auto push_vtx(int v, long weight) -> void
    {
        vv.push_back(v);
        total_wt += weight;
    }
    auto pop_vtx(long weight) -> void
    {
        total_wt -= weight;
        vv.pop_back();
    }
};

void add_edge(Graph *g, int v, int w);

void populate_bit_complement_nd(Graph & g);

Graph induced_subgraph(const Graph & g, const vector<int> vv);

#endif

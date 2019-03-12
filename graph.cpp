//#define _GNU_SOURCE

#include "graph.h"
#include "util.h"
#include "bitset.h"

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

Graph::Graph(int n)
        : n(n), numwords((n+BITS_PER_WORD-1)/BITS_PER_WORD), weight(n),
          bit_complement_nd(n, vector<unsigned long long>(numwords, ~0ull))
{
    if (0 != n % BITS_PER_WORD) {
        for (int i=0; i<n; i++) {
            bit_complement_nd[i].back() = (1ull << (n % BITS_PER_WORD)) - 1;
        }
    }
    for (int i=0; i<n; i++) {
        unset_bit(bit_complement_nd[i], i);
    }
}

auto Graph::resize(unsigned new_n) -> void
{
    n = new_n;
    numwords = ((n+BITS_PER_WORD-1)/BITS_PER_WORD);
    weight.resize(n);
    bit_complement_nd.resize(n);

    for (auto & bcn : bit_complement_nd) {
        bcn.assign(numwords, ~0ull);
    }

    if (0 != n % BITS_PER_WORD) {
        for (int i=0; i<n; i++) {
            bit_complement_nd[i].back() = (1ull << (n % BITS_PER_WORD)) - 1;
        }
    }
    for (int i=0; i<n; i++) {
        unset_bit(bit_complement_nd[i], i);
    }
}

Graph induced_subgraph(const Graph & g, const vector<int> vv) {
    struct Graph subg(vv.size());
    for (int i=0; i<subg.n; i++) {
        auto & row = g.bit_complement_nd[vv[i]];
        for (int j=i+1; j<subg.n; j++) {
            bool bit = test_bit(row, vv[j]);
            if (!bit) {
                unset_bit(subg.bit_complement_nd[i], j);
                unset_bit(subg.bit_complement_nd[j], i);
            }
        }
    }

    for (int i=0; i<subg.n; i++)
        subg.weight[i] = g.weight[vv[i]];
    return subg;
}

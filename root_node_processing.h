#ifndef ROOT_NODE_PROCESSING_H
#define ROOT_NODE_PROCESSING_H

#include "sparse_graph.h"

vector<int> initialise(const SparseGraph & g);

auto remove_vertices_with_closed_nd_wt_leq_incumbent(const SparseGraph & g,
        vector<int> & vv, long current_wt, long incumbent_wt, int num_threads) -> void;

auto reduce_and_reverse_cp_order(vector<int> & vv, const SparseGraph & g,
        long current_wt, long incumbent_wt) -> void;

#endif

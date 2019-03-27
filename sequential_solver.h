#ifndef COLOUR_ORDER_SOLVER_H
#define COLOUR_ORDER_SOLVER_H

#include "sparse_graph.h"
#include "graph.h"
#include "params.h"

#include <atomic>

auto sequential_mwc(const SparseGraph & g, const Params params, VtxList & incumbent, long & search_node_count) -> void;

#endif

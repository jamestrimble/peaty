#ifndef SPARSE_GRAPH_H
#define SPARSE_GRAPH_H

#include <limits.h>
#include <stdbool.h>
#include <stdlib.h>

#include <algorithm>
#include <iostream>
#include <thread>
#include <vector>

#include "graph.h"
#include "util.h"

using std::vector;
using Edge = std::pair<int, int>;

struct AdjListElement
{
    int v;
    int next;
};

struct ColourClass
{
    vector<int> vv;
    long weight;
};

struct SparseGraph
{
    unsigned n;  // number of vertices

    // an adjacency list for each vertex
    vector<vector<int>> adjlist;

    // vertex weights
    vector<long> weight;

    vector<bool> vertex_has_loop;

    SparseGraph(int n) : n(n), adjlist(n), weight(n, 1), vertex_has_loop(n)
    {
    }

    auto add_loop(int v) -> void
    {
        vertex_has_loop[v] = true;
    }

    auto add_edge(int v, int w) -> void
    {
        adjlist[v].push_back(w);
        adjlist[w].push_back(v);
    }

    auto has_edge(int v, int w) const -> bool
    {
        if (adjlist[w].size() < adjlist[v].size()) {
            // for speed, look in the smaller adjacency list
            std::swap(v, w);
        }
        auto & v_adj_lst = adjlist[v];
        return std::find(v_adj_lst.begin(), v_adj_lst.end(), w) != v_adj_lst.end();
    }

    auto vv_are_clique(const vector<int> & vv) const -> bool
    {
        for (unsigned i=0; i<vv.size(); i++) {
            int v = vv[i];
            auto & v_adj_lst = adjlist[v];
            for (unsigned j=i+1; j<vv.size(); j++) {
                int w = vv[j];
                if (std::find(v_adj_lst.begin(), v_adj_lst.end(), w) == v_adj_lst.end()) {
                    return false;
                }
            }
        }
        return true;
    }

    auto remove_edges_incident_to_loopy_vertices() -> void
    {
        for (unsigned i=0; i<n; i++) {
            if (vertex_has_loop[i]) {
                adjlist[i].clear();
            } else {
                auto & lst = adjlist[i];
                lst.erase(std::remove_if(
                            lst.begin(),
                            lst.end(),
                            [this](int v){return vertex_has_loop[v];}),
                        lst.end());
            }
        }
    }

    auto sort_adj_lists() -> void
    {
        for (auto & list : adjlist)
            std::sort(list.begin(), list.end());
    }

    auto complement_of_induced_subgraph(const vector<int> & vv) const -> Graph
    {
        vector<int> old_to_new_vtx(n, -1);
        Graph subgraph(vv.size());
        for (unsigned i=0; i<vv.size(); i++)
            old_to_new_vtx[vv[i]] = i;

        for (unsigned i=0; i<vv.size(); i++)
            for (unsigned j=0; j<vv.size(); j++)
                if (i != j)
                    subgraph.add_edge(i, j);

        for (int old_v : vv) {
            int new_v = old_to_new_vtx[old_v];
            for (int old_w : adjlist[old_v]) {
                if (old_w > old_v)
                    continue;
                int new_w = old_to_new_vtx[old_w];
                if (new_w != -1)
                    subgraph.remove_edge(new_v, new_w);
            }
        }
        for (unsigned i=0; i<vv.size(); i++)
            subgraph.weight[i] = weight[vv[i]];

        return subgraph;
    }

    template<typename T>
    auto induced_subgraph(const vector<int> & vv) const -> T
    {
        vector<int> old_to_new_vtx(n, -1);
        T subgraph(vv.size());
        induced_subgraph<T>(vv, old_to_new_vtx, subgraph);
        return subgraph;
    }

    template<typename T>
    auto induced_subgraph(const vector<int> & vv, vector<int> & old_to_new_vtx,
            T & subgraph) const -> void
    {
        for (unsigned i=0; i<vv.size(); i++)
            old_to_new_vtx[vv[i]] = i;

        for (int old_v : vv) {
            int new_v = old_to_new_vtx[old_v];
            for (int old_w : adjlist[old_v]) {
                if (old_w > old_v)
                    continue;
                int new_w = old_to_new_vtx[old_w];
                if (new_w != -1)
                    subgraph.add_edge(new_v, new_w);
            }
        }
        for (unsigned i=0; i<vv.size(); i++)
            subgraph.weight[i] = weight[vv[i]];

        for (unsigned i=0; i<vv.size(); i++)
            old_to_new_vtx[vv[i]] = -1;
    }

    template<typename T>
    auto induced_subgraph_reusing_graph(const vector<int> & vv, vector<int> & old_to_new_vtx,
            T & subgraph) const -> void
    {
        subgraph.resize(vv.size());
        induced_subgraph<T>(vv, old_to_new_vtx, subgraph);
    }

    auto resize(unsigned new_n) -> void
    {
        n = new_n;
        weight.resize(n);
        adjlist.resize(n);

        for (auto & lst : adjlist)
            lst.clear();
    }

    auto print_dimacs_format() -> void
    {
        long endpoint_count = 0;
        for (auto & lst : adjlist)
            endpoint_count += lst.size();
        std::cout << "p edge " << n << " " << (endpoint_count/2) << std::endl;
        for (unsigned i=0; i<n; i++) {
            std::cout << "n " << (i+1) << " " << weight[i] << std::endl;
        }
        for (unsigned i=0; i<n; i++) {
            for (unsigned j : adjlist[i]) {
                if (j > i)
                    break;
                std::cout << "e " << (i+1) << " " << (j+1) << std::endl;
            }
        }
    }
};

SparseGraph readSparseGraph();

SparseGraph readSparseGraphPaceFormat();

SparseGraph fastReadSparseGraphPaceFormat();

#endif

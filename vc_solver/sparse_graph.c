#include "sparse_graph.h"
#include "util.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <ios>
#include <iostream>
#include <sstream>

auto deduplicate_and_add_edges(SparseGraph & g, vector<Edge> & edges)
{
    std::sort(edges.begin(), edges.end());
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
    printf("c %d edges after de-duplication\n", int(edges.size()));

    for (Edge e : edges)
        g.add_edge(e.first, e.second);

    g.sort_adj_lists();
}

SparseGraph readSparseGraph(char* filename) {
    int medges = 0;
    int edges_read = 0;

    struct SparseGraph g(0);

    vector<Edge> edges;

    {
        std::string line;
        while (std::getline(std::cin, line)) {
            std::stringstream line_stream(line);
            std::string token;
            if (!(line_stream >> token))
                continue;
            if (token == "e") {
                int v, w;
                if (!(line_stream >> v) || !(line_stream >> w))
                    fail("Error reading a line beginning with e.\n");
                if (v < w)
                    edges.push_back({v-1, w-1});
                else if (v > w)
                    edges.push_back({w-1, v-1});
                edges_read++;
            } else if (token == "p") {
                long n, m;
                if (!(line_stream >> token) || token != "edge" || !(line_stream >> n) || !(line_stream >> m))
                    fail("Error reading a line beginning with p.\n");
                printf("c %ld vertices\n", n);
                printf("c %ld edges\n", m);
                if (n > INT_MAX)
                    fail("Too many vertices.\n");
                if (m > INT_MAX)
                    fail("Too many edges.\n");
                medges = m;
                g = SparseGraph(n);
            } else if (token == "n") {
                int v;
                long wt;
                if (!(line_stream >> v) || !(line_stream >> wt))
                    fail("Error reading a line beginning with n.\n");
                g.weight[v-1] = wt;
            }
        }
    }

    if (medges>0 && int(edges.size()) != medges)
        fail("Unexpected number of edges.");

    deduplicate_and_add_edges(g, edges);

    return g;
}

SparseGraph readSparseGraphPaceFormat(char* filename) {
    int medges = 0;
    int edges_read = 0;

    struct SparseGraph g(0);

    vector<Edge> edges;

    {
        std::string line;
        while (std::getline(std::cin, line)) {
            std::stringstream line_stream(line);
            std::string token;
            if (!(line_stream >> token))
                continue;
            if (token[0] >= '0' && token[0] <= '9') {
                int v=-1, w;
                try {
                    v = std::stoi(token);
                } catch (std::invalid_argument) {
                    fail("Invalid edge");
                }
                if (!(line_stream >> w))
                    fail("Error reading a line beginning with e.\n");
                if (v < w)
                    edges.push_back({v-1, w-1});
                else if (v > w)
                    edges.push_back({w-1, v-1});
                edges_read++;
            } else if (token == "p") {
                long n, m;
                if (!(line_stream >> token) || token != "td" || !(line_stream >> n) || !(line_stream >> m))
                    fail("Error reading a line beginning with p.\n");
                printf("c %ld vertices\n", n);
                printf("c %ld edges\n", m);
                if (n > INT_MAX)
                    fail("Too many vertices.\n");
                if (m > INT_MAX)
                    fail("Too many edges.\n");
                medges = m;
                g = SparseGraph(n);
            }
        }
    }

    if (medges>0 && int(edges.size()) != medges)
        fail("Unexpected number of edges.");

    deduplicate_and_add_edges(g, edges);

    return g;
}

#include "sparse_graph.h"
#include "util.h"

#include <algorithm>
#include <cctype>
//#include <fstream>
#include <ios>
#include <iostream>
#include <sstream>

auto deduplicate_and_add_edges(SparseGraph & g, vector<Edge> & edges)
{
    std::sort(edges.begin(), edges.end());
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
    printf("c %d edges after de-duplication\n", int(edges.size()));

    for (Edge e : edges) {
        int v = e.first;
        int w = e.second;
        if (v == w) {
            g.add_loop(v);
        } else {
            g.add_edge(v, w);
        }
    }

    g.sort_adj_lists();
}

SparseGraph readSparseGraph() {
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

SparseGraph readSparseGraphPaceFormat() {
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
                else
                    edges.push_back({w-1, v-1});
                edges_read++;
                if (edges_read == medges)
                    break;
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

struct FastFileReader
{
private:
    std::istream & file;
    char c;

    void advance()
    {
        if (!(file >> c))
            c = 0;
    }

public:
    FastFileReader(std::istream & file) : file(file)
    {
        file.unsetf(std::ios_base::skipws);
        advance();
    }

    int get_digit()
    {
        char ch = c;
        advance();
        return ch - '0';
    }

    long get_long()
    {
        long retval = 0;
        while (is_digit()) {
            retval *= 10;
            retval += get_digit();
        }
        return retval;
    }

    bool is_digit()
    {
        return c >= '0' && c <= '9';
    }

    bool finished()
    {
        return c == 0;
    }

    char peek_char()
    {
        return c;
    }

    char get_char()
    {
        char ch = c;
        advance();
        return ch;
    }

    void skip_ws()
    {
        while (std::isspace(c))
            advance();
    }

    void skip_ws_alpha()
    {
        while (std::isspace(c) || std::isalpha(c))
            advance();
    }

    void skip_to_end_of_line()
    {
        while (c != 10 && c != 13)
            advance();
        skip_ws();
    }
};

SparseGraph fastReadSparseGraphPaceFormat() {
    int medges = 0;
    int edges_read = 0;

    struct SparseGraph g(0);

    vector<Edge> edges;

    {
        FastFileReader ffr(std::cin);

        std::string line;

        while (!ffr.finished()) {
            char ch = ffr.peek_char();
            if (ch == 'c') {
                ffr.skip_to_end_of_line();
            } else if (ffr.is_digit()) {
                int v = ffr.get_long();
                ffr.skip_ws();
                if (!ffr.is_digit())
                    fail("Error reading an edge line.\n");
                int w = ffr.get_long();
                ffr.skip_ws();
                if (v < w)
                    edges.push_back({v-1, w-1});
                else
                    edges.push_back({w-1, v-1});
                edges_read++;
            } else if (ch == 'p') {
                ffr.skip_ws_alpha();
                if (!ffr.is_digit())
                    fail("Error reading a line beginning with p.\n");
                long n = ffr.get_long();
                ffr.skip_ws();
                if (!ffr.is_digit())
                    fail("Error reading a line beginning with p.\n");
                long m = ffr.get_long();
                ffr.skip_ws();
                printf("c %ld vertices\n", n);
                printf("c %ld edges\n", m);
                if (n > INT_MAX)
                    fail("Too many vertices.\n");
                if (m > INT_MAX)
                    fail("Too many edges.\n");
                medges = m;
                g = SparseGraph(n);
            } else {
                fail("Unexpected character in input.");
            }
        }
    }

    if (medges>0 && int(edges.size()) != medges)
        fail("Unexpected number of edges.");

    deduplicate_and_add_edges(g, edges);

    return g;
}

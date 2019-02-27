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

struct FastFileReader
{
private:
    std::ifstream & file;
    char c;

    void advance()
    {
        if (!(file >> c))
            c = 0;
    }

public:
    FastFileReader(std::ifstream & file) : file(file)
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

SparseGraph fastReadSparseGraph(char* filename) {
    int medges = 0;
    int edges_read = 0;

    struct SparseGraph g(0);

    vector<Edge> edges;

    {
        std::ifstream file(filename);
        if (!file.good())
            fail("Cannot open file");

        FastFileReader ffr(file);

        std::string line;

        while (!ffr.finished()) {
            char ch = ffr.get_char();
            if (ch == 'c') {
                ffr.skip_to_end_of_line();
            } else if (ch == 'e') {
                ffr.skip_ws();
                if (!ffr.is_digit())
                    fail("Error reading a line beginning with e.\n");
                int v = ffr.get_long();
                ffr.skip_ws();
                if (!ffr.is_digit())
                    fail("Error reading a line beginning with e.\n");
                int w = ffr.get_long();
                ffr.skip_ws();
                if (v < w)
                    edges.push_back({v-1, w-1});
                else if (v > w)
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
            } else if (ch == 'n' || ch == 'v') {
                ffr.skip_ws();
                if (!ffr.is_digit())
                    fail("Error reading a line beginning with n or v.\n");
                int v = ffr.get_long();
                ffr.skip_ws();
                if (!ffr.is_digit())
                    fail("Error reading a line beginning with n or v.\n");
                long wt= ffr.get_long();
                ffr.skip_ws();
                g.weight[v-1] = wt;
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

SparseGraph readSparseGraphEdgesFormat(char* filename, unsigned num_lines_to_skip) {
    struct SparseGraph g(0);

    vector<Edge> edges;

    {
        std::ifstream file(filename);
        if (!file.good())
            fail("Cannot open file");

        long v, w;
        std::string line;
        unsigned long long num_lines = 0;
        while (std::getline(file, line)) {
            if (line.substr(0, 1) == "%")
                continue;
            ++num_lines;
            if (num_lines <= num_lines_to_skip)
                continue;
            std::stringstream line_stream(line);
            if (!(line_stream >> v) || !(line_stream >> w))
                fail("Error reading an edge.\n");

            if (v > INT_MAX || w > INT_MAX)
                fail("A vertex number is too large.\n");
            if (v < w)
                edges.push_back({v-1, w-1});
            else if (v > w)
                edges.push_back({w-1, v-1});
        }
    }

    if (edges.size() > unsigned(INT_MAX))
        fail("Too many edges.\n");

    if (edges.size()) {
        int min_v = INT_MAX;
        int max_v = 0;
        for (Edge e : edges) {
            if (e.second > max_v)
                max_v = e.second;
            if (e.first < min_v)
                min_v = e.first;
        }
        if (min_v == -1) {
            ++min_v;
            ++max_v;
            for (Edge & e : edges) {
                ++e.first;
                ++e.second;
            }
        }
        int n = max_v + 1;
        g = SparseGraph(n);
        for (int i=0; i<n; i++)
            g.weight[i] = 1 + ((i+1) % 200);

        deduplicate_and_add_edges(g, edges);
    }

    return g;
}

SparseGraph readSparseGraphMtxFormat(char* filename) {
    int medges = 0;

    struct SparseGraph g(0);

    vector<Edge> edges;

    int num_edges_read = 0;
    {
        std::ifstream file(filename);
        if (!file.good())
            fail("Cannot open file");

        std::string line;
        bool read_header = false;
        while (std::getline(file, line)) {
            std::stringstream line_stream(line);
            std::string token;
            if (!(line_stream >> token))
                continue;
            if (token.substr(0, 1) == "%")
                continue;
            long n1, n2, m;
            n1 = atol(token.c_str());
            if (!(line_stream >> n2) || !(line_stream >> m))
                fail("Error reading line that states size of matrix.\n");
            printf("c n1 = %ld\n", n1);
            printf("c n2 = %ld\n", n2);
            printf("c m = %ld edges\n", m);
            long n = std::max(n1, n2);
            if (n > INT_MAX)
                fail("Too many vertices.\n");
            if (m > INT_MAX)
                fail("Too many edges.\n");
            medges = m;
            edges.reserve(m);
            g = SparseGraph(n);
            for (int i=0; i<n; i++)
                g.weight[i] = 1 + ((i+1) % 200);
            read_header = true;
            break;
        }
        if (!read_header)
            fail("Failed to read header.\n");

        int v, w;
        while (std::getline(file, line)) {
            if (line.substr(0, 1) == "%")
                continue;
            std::stringstream line_stream(line);
            if (!(line_stream >> v) || !(line_stream >> w))
                fail("Error reading an edge.\n");

            if (v < w)
                edges.push_back({v-1, w-1});
            else if (v > w)
                edges.push_back({w-1, v-1});
            ++num_edges_read;
        }
    }

    // num_edges_read is used instead of edges.size() because it includes (ignored) loops
    if (medges>0 && num_edges_read != medges)
        fail("Unexpected number of edges.");

    deduplicate_and_add_edges(g, edges);

    return g;
}

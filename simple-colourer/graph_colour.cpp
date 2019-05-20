#define _POSIX_SOURCE

#include "../graph_colour_solver.h"

#include <iostream>
#include <string>
#include <vector>

#include <argp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void fail(const std::string & msg) {
    std::cerr << msg << std::endl;
    exit(1);
}


static char doc[] = "Colour a graph read from standard input.  The expected format is the DIMACS clique format.";
static char args_doc[] = "FILENAME";
static struct argp_option options[] = {
    {"fractional-level", 'f', "LEVEL", 0, "1 for colouring, 2 for two colours per vertex, etc"},
    {"time-limit", 'l', "LIMIT", 0, "Time limit in milliseconds"},
    { 0 }
};

static struct {
    int fractional_level;
    int time_limit;
    int arg_num;
} arguments;

void set_default_arguments() {
    arguments.fractional_level = 1;
    arguments.time_limit = 0;
    arguments.arg_num = 0;
}

static error_t parse_opt (int key, char *arg, struct argp_state *state) {
    switch (key) {
        case 'f':
            arguments.fractional_level = atoi(arg);
            break;
        case 'l':
            arguments.time_limit = atoi(arg);
            break;
        case ARGP_KEY_ARG:
            argp_usage(state);
            break;
        case ARGP_KEY_END:
            break;
        default: return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

static struct argp argp = { options, parse_opt, args_doc, doc };

/******************************************************************************/

struct ColouringGraph readColouringGraph()
{
    char* line = NULL;
    size_t nchar = 0;

    int nvertices = 0;
    int medges = 0;
    int v, w;
    int edges_read = 0;

    struct ColouringGraph g(0);

    while (getline(&line, &nchar, stdin) != -1) {
        if (nchar > 0) {
            switch (line[0]) {
            case 'p':
                if (sscanf(line, "p edge %d %d", &nvertices, &medges)!=2)
                    fail("Error reading a line beginning with p.\n");
                printf("%d vertices\n", nvertices);
                printf("%d edges\n", medges);
                g = ColouringGraph(nvertices);
                break;
            case 'e':
                if (sscanf(line, "e %d %d", &v, &w)!=2)
                    fail("Error reading a line beginning with e.\n");
                v -= 1; // since input format is 1-based
                w -= 1; // ditto
                g.add_edge(v, w);
                edges_read++;
                break;
            }
        }
    }

    if (medges>0 && edges_read != medges) fail("Unexpected number of edges.");

    return g;
}

int main(int argc, char** argv)
{
    set_default_arguments();
    argp_parse(&argp, argc, argv, 0, 0, 0);

    struct ColouringGraph g = readColouringGraph();
    struct ColouringGraph complement_g(g.n);
    for (int i=0; i<g.n; i++) {
        for (int j=i+1; j<g.n; j++) {
            if (!g.adj_matrix[i][j]) {
                complement_g.add_edge(i, j);
            }
        }
    }

    std::atomic_bool terminate_early(false);
    int colouring_number = find_colouring_number(complement_g, arguments.fractional_level, terminate_early);

    printf("%d-fold colouring number is %d\n", arguments.fractional_level, colouring_number);
}

#define _POSIX_SOURCE

#include "c_program_timing.h"
#include "graph_colour_solver.h"

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


static char doc[] = "Colour a graph.  The expected file format is the DIMACS clique format.";
static char args_doc[] = "FILENAME";
static struct argp_option options[] = {
    {"quiet", 'q', 0, 0, "Quiet output"},
    {"rng-seed", 'r', "SEED", 0, "seed for RNG to add a hint of randomness to vertex order"},
    {"fractional-level", 'f', "LEVEL", 0, "1 for colouring, 2 for two colours per vertex, etc"},
    {"time-limit", 'l', "LIMIT", 0, "Time limit in milliseconds"},
    { 0 }
};

static struct {
    bool quiet;
    int rng_seed;
    int fractional_level;
    int time_limit;
    char *filename;
    int arg_num;
} arguments;

void set_default_arguments() {
    arguments.quiet = false;
    arguments.rng_seed = 1;
    arguments.fractional_level = 1;
    arguments.time_limit = 0;
    arguments.filename = NULL;
    arguments.arg_num = 0;
}

static error_t parse_opt (int key, char *arg, struct argp_state *state) {
    switch (key) {
        case 'q':
            arguments.quiet = true;
            break;
        case 'r':
            arguments.rng_seed = atoi(arg);
            break;
        case 'f':
            arguments.fractional_level = atoi(arg);
            break;
        case 'l':
            arguments.time_limit = atoi(arg);
            break;
        case ARGP_KEY_ARG:
            if (arguments.arg_num >= 1)
                argp_usage(state);
            arguments.filename = arg;
            arguments.arg_num++;
            break;
        case ARGP_KEY_END:
            if (arguments.arg_num == 0)
                argp_usage(state);
            break;
        default: return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

static struct argp argp = { options, parse_opt, args_doc, doc };

/******************************************************************************/

struct Graph readGraph(char* filename, int fractional_level)
{
    char* line = NULL;
    size_t nchar = 0;

    int nvertices = 0;
    int medges = 0;
    int v, w;
    int edges_read = 0;

    struct Graph g(0);

    while (getline(&line, &nchar, stdin) != -1) {
        if (nchar > 0) {
            switch (line[0]) {
            case 'p':
                if (sscanf(line, "p edge %d %d", &nvertices, &medges)!=2)
                    fail("Error reading a line beginning with p.\n");
                printf("%d vertices\n", nvertices);
                printf("%d edges\n", medges);
                g = Graph(nvertices);
                break;
            case 'e':
                if (sscanf(line, "e %d %d", &v, &w)!=2)
                    fail("Error reading a line beginning with e.\n");
                v -= 1; // since input format is 1-based
                w -= 1; // ditto
                add_edge(g, v, w);
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

    struct Graph input_g_unsorted = readGraph(arguments.filename, arguments.fractional_level);

    std::vector<int> shuffled_ints;
    for (int i=0; i<input_g_unsorted.n; i++)
        shuffled_ints.push_back(i);
    srand(arguments.rng_seed);
    for (int i=input_g_unsorted.n-1; i>=1; i--) {
        int r = rand() % (i+1);
        std::swap(shuffled_ints[i], shuffled_ints[r]);
    }

    std::vector<int> vv;

    std::vector<int> degree(input_g_unsorted.n);
    for (int v=0; v<input_g_unsorted.n; v++) {
        for (int w=v+1; w<input_g_unsorted.n; w++) {
            if (input_g_unsorted.adj_matrix[v][w]) {
                ++degree[v];
                ++degree[w];
            }
        }
    }

    for (int deg=0; deg<input_g_unsorted.n; deg++)
        for (int v : shuffled_ints)
            if (degree[v] == deg)
                vv.push_back(v);
    struct Graph input_g = induced_subgraph(&input_g_unsorted, vv);

    set_start_time();
    set_time_limit_ms(arguments.time_limit);

    for (int num_colours=0; ; num_colours++) {
        struct Solution clq;
        init_Solution(&clq, input_g.n);

        long expand_call_count = 0;
        make_adjacency_lists(&input_g);
        solve(&input_g, &expand_call_count, arguments.quiet, &clq, num_colours, arguments.fractional_level);
        long elapsed_msec = get_elapsed_time_msec();
        if (is_timeout_flag_set()) {
            elapsed_msec = arguments.time_limit * 1000;
            printf("%s %d %ld %s\n", arguments.filename, num_colours, elapsed_msec, "TIMEOUT");
            break;
        }

//        if (clq.size == input_g->n * arguments.fractional_level) {
//            printf("Solution");
//            struct Solution solution_for_unsorted_graph;
//            init_Solution(&solution_for_unsorted_graph, input_g->n);
//            solution_for_unsorted_graph.size = clq.size;
//            for (int i=0; i<input_g->n; i++) {
//                solution_for_unsorted_graph.vtx_colour[vv[i]] = clq.vtx_colour[i];
//            }
//            if (!is_solution_valid(input_g_unsorted, &solution_for_unsorted_graph, num_colours))
//                fail("The solution that was found is not a colouring with the correct number of colours!!!");
//            for (int i=0; i<input_g->n; i++) {
//                printf(" %d", solution_for_unsorted_graph.vtx_colour[i]);
//            }
//            destroy_Solution(&solution_for_unsorted_graph);
//            printf("\n");
//        }

        printf("%d %ld %s\n", num_colours, expand_call_count, clq.size == input_g.n * arguments.fractional_level ? "SATISFIABLE" : "UNSAT");

        destroy_Solution(&clq);

        if (clq.size == input_g.n * arguments.fractional_level)
            break;
    }
}

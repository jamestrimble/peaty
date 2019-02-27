//#define _GNU_SOURCE
#define _POSIX_SOURCE

#include "sparse_graph.h"
#include "util.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <algorithm>

int main(int argc, char** argv) {
    SparseGraph g = readSparseGraphEdgesFormat(argv[1], atoi(argv[2]));
    g.print_dimacs_format();
}

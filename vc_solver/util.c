#include <stdio.h>
#include <stdlib.h>

void fail(const char* msg) {
    fprintf(stderr, "%s\n", msg);
    exit(1);
}



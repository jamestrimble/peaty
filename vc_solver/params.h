#ifndef PARAMS_H
#define PARAMS_H

struct Params
{
    int colouring_variant;
    int max_sat_level;
    int algorithm_num;
    int num_threads;
    bool quiet;
    bool unweighted_sort;

    Params(int colouring_variant, int max_sat_level, int algorithm_num, int num_threads,
            bool quiet, int unweighted_sort) :
            colouring_variant(colouring_variant),
            max_sat_level(max_sat_level),
            algorithm_num(algorithm_num),
            num_threads(num_threads),
            quiet(quiet),
            unweighted_sort(unweighted_sort)
    {}
};

#endif

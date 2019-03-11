#ifndef PARAMS_H
#define PARAMS_H

#include<atomic>
#include<chrono>

struct Params
{
    int colouring_variant;
    int max_sat_level;
    int algorithm_num;
    int num_threads;
    bool quiet;
    bool unweighted_sort;
    int ind_set_upper_bound;

    Params(int colouring_variant, int max_sat_level, int algorithm_num, int num_threads,
            bool quiet, int unweighted_sort, int ind_set_upper_bound) :
            colouring_variant(colouring_variant),
            max_sat_level(max_sat_level),
            algorithm_num(algorithm_num),
            num_threads(num_threads),
            quiet(quiet),
            unweighted_sort(unweighted_sort),
            ind_set_upper_bound(ind_set_upper_bound)
    {}
};

#endif

#define _POSIX_C_SOURCE 199309L

#include <stdbool.h>
#include <time.h>
#include "c_program_timing.h"

struct timespec start_time;

long time_limit_sec = 0;

bool timed_out = false;

// Based on https://gist.github.com/diabloneo/9619917
struct timespec timespec_diff(struct timespec start,
                              struct timespec stop)
{
    if (stop.tv_nsec < start.tv_nsec) {
        stop.tv_nsec += 1000000000;
        stop.tv_sec--;
    }
    struct timespec result;
    result.tv_sec = stop.tv_sec - start.tv_sec;
    result.tv_nsec = stop.tv_nsec - start.tv_nsec;
    return result;
}

long timespec_to_msec(struct timespec ts)
{
    return ts.tv_sec*1000 + ts.tv_nsec/1000000;
}

void set_time_limit_sec(long time_limit)
{
    time_limit_sec = time_limit;
}

void check_for_timeout()
{
    timed_out = time_limit_sec && get_elapsed_timespec().tv_sec >= time_limit_sec;
}

bool is_timeout_flag_set()
{
    return timed_out;
}

void set_start_time()
{
    clock_gettime(CLOCK_MONOTONIC, &start_time);
}

struct timespec get_elapsed_timespec()
{
    struct timespec end_time;
    clock_gettime(CLOCK_MONOTONIC, &end_time);
    return timespec_diff(start_time, end_time);
}

long get_elapsed_time_msec()
{
    return timespec_to_msec(get_elapsed_timespec());
}


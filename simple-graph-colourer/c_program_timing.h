#include <stdbool.h>

void set_start_time();

struct timespec get_elapsed_timespec();

long get_elapsed_time_msec();

void set_time_limit_sec(long time_limit);

void check_for_timeout();

bool is_timeout_flag_set();

#define _POSIX_SOURCE

#include "graph.h"
#include "sparse_graph.h"
#include "util.h"
#include "sequential_solver.h"
#include "params.h"

#include <argp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <functional>
#include <mutex>
#include <thread>

using std::atomic;
using std::condition_variable;
using std::cv_status;
using std::function;
using std::mutex;
using std::thread;
using std::unique_lock;

using std::chrono::seconds;
using std::chrono::steady_clock;
using std::chrono::system_clock;
using std::chrono::duration_cast;
using std::chrono::milliseconds;

static char doc[] = "Find a maximum clique in a graph in DIMACS format";
static char args_doc[] = "FILENAME";
static struct argp_option options[] = {
    {"quiet", 'q', 0, 0, "Quiet output"},
    {"unweighted-sort", 'u', 0, 0, "Unweighted ordering (only applies to certain algorithms)"},
    {"colouring-variant", 'c', "VARIANT", 0, "For algorithms 0 and 5, which type of colouring?"},
    {"ind-set-upper-bound", 'i', "BOUND", 0, "An upper bound on the size of the max independent set"},
    {"time-limit", 'l', "LIMIT", 0, "Time limit in seconds"},
    {"algorithm", 'a', "NUMBER", 0, "Algorithm number"},
    {"max-sat-level", 'm', "LEVEL", 0, "Level of MAXSAT reasoning; default=2"},
    {"num-threads", 't', "NUMBER", 0, "Number of threads (for parallel algorithms only)"},
    {"file-format", 'f', "FORMAT", 0, "File format (DIMACS, MTX or EDGES)"},
    { 0 }
};

enum class FileFormat
{
    Dimacs,
    Pace
};

static struct {
    bool quiet;
    bool unweighted_sort;
    int colouring_variant;
    int ind_set_upper_bound;
    int time_limit;
    int algorithm_num;
    int max_sat_level;
    int num_threads;
    FileFormat file_format;
    char *filename;
    int arg_num;
} arguments;

void set_default_arguments()
{
    arguments.file_format = FileFormat::Dimacs;
    arguments.num_threads = 1;
}

static error_t parse_opt (int key, char *arg, struct argp_state *state)
{
    switch (key) {
        case 'q':
            arguments.quiet = true;
            break;
        case 'u':
            arguments.unweighted_sort = true;
            break;
        case 'c':
            arguments.colouring_variant = atoi(arg);
            break;
        case 'i':
            arguments.ind_set_upper_bound = atoi(arg);
            break;
        case 'l':
            arguments.time_limit = atoi(arg);
            break;
        case 'a':
            arguments.algorithm_num = atoi(arg);
            break;
        case 'm':
            arguments.max_sat_level = atoi(arg);
            break;
        case 't':
            arguments.num_threads = atoi(arg);
            break;
        case 'f':
            if (!strcmp(arg, "PACE") || !strcmp(arg, "pace"))
                arguments.file_format = FileFormat::Pace;
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

// Checks if a set of vertices induces a clique
bool check_indset(const SparseGraph & g, const VtxList & clq) {
    long total_wt = 0;
    for (unsigned i=0; i<clq.vv.size(); i++)
        total_wt += g.weight[clq.vv[i]];
    if (total_wt != clq.total_wt)
        return false;

    for (unsigned i=0; i<clq.vv.size(); i++) {
        int v = clq.vv[i];
        for (unsigned j=i+1; j<clq.vv.size(); j++) {
            int w = clq.vv[j];
            if (v == w) {
                return false;
            }
            if (std::find(g.adjlist[v].begin(), g.adjlist[v].end(), w) != g.adjlist[v].end()) {
                return false;
            }
        }
    }
    return true;
}


// Credit for timing and timeout code: Ciaran's code from our paper on restarting SIP

/* Helper: return a function that runs the specified algorithm, dealing
 * with timing information and timeouts. */
template <typename Result_, typename Params_, typename Data_>
auto run_this_wrapped(const function<Result_ (const Data_ &, const Params_ &)> & func)
    -> function<Result_ (const Data_ &, Params_ &, bool &, int)>
{
    return [func] (const Data_ & data, Params_ & params, bool & aborted, int timeout) -> Result_ {
        /* For a timeout, we use a thread and a timed CV. We also wake the
         * CV up if we're done, so the timeout thread can terminate. */
        thread timeout_thread;
        mutex timeout_mutex;
        condition_variable timeout_cv;
        atomic<bool> abort;
        abort.store(false);
        params.abort = &abort;
        if (0 != timeout) {
            timeout_thread = thread([&] {
                    auto abort_time = steady_clock::now() + seconds(timeout);
                    {
                        /* Sleep until either we've reached the time limit,
                         * or we've finished all the work. */
                        unique_lock<mutex> guard(timeout_mutex);
                        while (! abort.load()) {
                            if (cv_status::timeout == timeout_cv.wait_until(guard, abort_time)) {
                                /* We've woken up, and it's due to a timeout. */
                                aborted = true;
                                break;
                            }
                        }
                    }
                    abort.store(true);
                    });
        }

        /* Start the clock */
        params.start_time = steady_clock::now();
        auto result = func(data, params);

        /* Clean up the timeout thread */
        if (timeout_thread.joinable()) {
            {
                unique_lock<mutex> guard(timeout_mutex);
                abort.store(true);
                timeout_cv.notify_all();
            }
            timeout_thread.join();
        }

        return result;
    };
}

/* Helper: return a function that runs the specified algorithm, dealing
 * with timing information and timeouts. */
template <typename Result_, typename Params_, typename Data_>
auto run_this(Result_ func(const Data_ &, const Params_ &)) -> function<Result_ (const Data_ &, Params_ &, bool &, int)>
{
    return run_this_wrapped(function<Result_ (const Data_ &, const Params_ &)>(func));
}

struct Result
{
    VtxList clq;
    long search_node_count;
    Result(const SparseGraph & g) : clq(g.n), search_node_count(0) {}
};

auto mwc(const SparseGraph & g, const Params & params) -> Result
{
    Result result(g);
    sequential_mwc(g, params, result.clq, result.search_node_count);
    return result;
}

int main(int argc, char** argv) {
    set_default_arguments();
    argp_parse(&argp, argc, argv, 0, 0, 0);

    if (arguments.algorithm_num != 5)
        arguments.num_threads = 1;

    const SparseGraph g =
            arguments.file_format==FileFormat::Pace ? readSparseGraphPaceFormat(arguments.filename) :
                                                      readSparseGraph(arguments.filename);

    Params params {arguments.colouring_variant, arguments.max_sat_level, arguments.algorithm_num,
            arguments.num_threads, arguments.quiet, arguments.unweighted_sort, arguments.ind_set_upper_bound};

    bool aborted = false;
    Result result = run_this(mwc)(g, params, aborted, arguments.time_limit);

    auto elapsed_msec = duration_cast<milliseconds>(steady_clock::now() - params.start_time).count();

    if (aborted) {
        printf("TIMEOUT\n");
        elapsed_msec = arguments.time_limit * 1000;
    }

    // sort vertices in clique by index
    std::sort(result.clq.vv.begin(), result.clq.vv.end());

    std::vector<bool> in_clq(g.n);
    for (unsigned i=0; i<result.clq.vv.size(); i++)
        in_clq[result.clq.vv[i]] = true;
    printf("VertexCover ");
    for (unsigned i=0; i<g.n; i++)
        if (!in_clq[i])
            printf("%d ", i);
    printf("\n");

    printf("Stats: status filename program algorithm_number max_sat_level num_threads size weight time_ms nodes\n");
    std::string fname(arguments.filename);
    auto last_slash = fname.find_last_of("/");
    if (last_slash < fname.size())
        fname = fname.substr(last_slash + 1, fname.size());
    std::cout <<
            (aborted ? "TIMEOUT" : "COMPLETED") << " " <<
            fname << " " <<
            argv[0] << " " <<
            arguments.algorithm_num << " " <<
            arguments.max_sat_level << " " <<
            arguments.num_threads << " " <<
            result.clq.vv.size() << " " <<
            result.clq.total_wt <<  " " <<
            elapsed_msec << " " <<
            result.search_node_count << std::endl;

    if (!check_indset(g, result.clq))
        fail("*** Error: the set of vertices found do not induce a clique of the expected weight\n");
}

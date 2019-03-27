#define _POSIX_SOURCE

#include "graph_colour_solver.h"

#include <limits.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

////////////////////////////////////////////////////////////////////////////////
//                                GRAPH STUFF                                 //
////////////////////////////////////////////////////////////////////////////////

void ColouringGraph::make_adjacency_lists()
{
    for (int i=0; i<n; i++)
        for (int j=0; j<n; j++)
            if (adj_matrix[i][j])
                adjlist[i].push_back(j);
}

ColouringGraph ColouringGraph::induced_subgraph(std::vector<int> & vv) const
{
    ColouringGraph subg(vv.size());
    for (int i=0; i<subg.n; i++)
        for (int j=0; j<i; j++)
            if (adj_matrix[vv[i]][vv[j]])
                subg.add_edge(i, j);

    return subg;
}

////////////////////////////////////////////////////////////////////////////////
//                                BITSET STUFF                                //
////////////////////////////////////////////////////////////////////////////////

#define BYTES_PER_WORD sizeof(unsigned long long)
#define BITS_PER_WORD (CHAR_BIT * BYTES_PER_WORD)

static bool test_bit(unsigned long long *bitset, int bit)
{
    return 0 != (bitset[bit/BITS_PER_WORD] & (1ull << (bit%BITS_PER_WORD)));
}

static void set_first_n_bits(unsigned long long *bitset, int n)
{
    for (int bit=0; bit<n; bit++)
        bitset[bit/BITS_PER_WORD] |= (1ull << (bit%BITS_PER_WORD));
}

static void unset_bit(unsigned long long *bitset, int bit)
{
    bitset[bit/BITS_PER_WORD] &= ~(1ull << (bit%BITS_PER_WORD));
}

static int bitset_popcount(unsigned long long *bitset, int num_words)
{
    int count = 0;
    for (int i=num_words-1; i>=0; i--)
        count += __builtin_popcountll(bitset[i]);
    return count;
}

static int bitset_intersection_popcount(unsigned long long *bitset1, unsigned long long *bitset2, int num_words)
{
    int count = 0;
    for (int i=num_words-1; i>=0; i--)
        count += __builtin_popcountll(bitset1[i] & bitset2[i]);
    return count;
}

static void bitset_intersect_with(unsigned long long *bitset1, unsigned long long *bitset2, int num_words)
{
    for (int i=num_words-1; i>=0; i--)
        bitset1[i] &= bitset2[i];
}

static bool bitset_empty(unsigned long long *bitset, int num_words)
{
    for (int i=num_words-1; i>=0; i--)
        if (bitset[i])
            return false;
    return true;
}

static void clear_bitset(unsigned long long *bitset, int num_words)
{
    for (int i=num_words-1; i>=0; i--)
        bitset[i] = 0ull;
}

static int first_set_bit(unsigned long long *bitset,
                         int num_words)
{
    for (int i=0; i<num_words; i++)
        if (bitset[i] != 0)
            return i*BITS_PER_WORD + __builtin_ctzll(bitset[i]);
    return -1;
}

static void copy_bitset(unsigned long long *src,
                        unsigned long long *dest,
                        int num_words)
{
    for (int i=0; i<num_words; i++)
        dest[i] = src[i];
}

////////////////////////////////////////////////////////////////////////////////
//                               SOLUTION STUFF                               //
////////////////////////////////////////////////////////////////////////////////

void solution_colour_vtx(struct Solution *solution, int v, int colour,
        unsigned long long *available_classes_bitset, int *num_colours_assigned_to_vertex, int domain_num_words, int f)
{
    ++solution->size;
    solution->vtx_colour[v][num_colours_assigned_to_vertex[v]] = colour;
    ++num_colours_assigned_to_vertex[v];
    unset_bit(available_classes_bitset + v * domain_num_words, colour);
    if (num_colours_assigned_to_vertex[v] == f)
        clear_bitset(available_classes_bitset + v * domain_num_words, domain_num_words);
}

void solution_pop_vtx(struct Solution *solution)
{
    solution->size--;
}

void solution_resize(struct Solution *solution, int size)
{
    solution->size = size;
}

void copy_Solution(struct Solution *src, struct Solution *dest)
{
    dest->size = src->size;
    dest->vtx_colour = src->vtx_colour;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// Precondition: at least one vertex remains that can be branched on.
int choose_branching_vertex(ColouringGraph *g, unsigned long long *available_classes_bitset,
        int domain_num_words)
{
    int best_available_class_count = INT_MAX;
    std::vector<int> vertices_with_best_available_class_count;
    for (int i=0; i<g->n; i++) {
        if (bitset_empty(available_classes_bitset + i * domain_num_words, domain_num_words))
            continue;
        int available_class_count = bitset_popcount(available_classes_bitset + i * domain_num_words, domain_num_words);
        if (available_class_count < best_available_class_count) {
            best_available_class_count = available_class_count;
            vertices_with_best_available_class_count.clear();
        }
        if (available_class_count == best_available_class_count)
            vertices_with_best_available_class_count.push_back(i);
    }

    std::vector<int> scores(vertices_with_best_available_class_count.size());

    for (unsigned i=0; i<vertices_with_best_available_class_count.size(); i++) {
        int v = vertices_with_best_available_class_count[i];
        for (unsigned j=0; j<i; j++) {
            int w = vertices_with_best_available_class_count[j];
            if (!g->adj_matrix[v][w]) {
                int pc = bitset_intersection_popcount(
                        available_classes_bitset + v * domain_num_words, available_classes_bitset + w * domain_num_words,
                        domain_num_words);
                scores[i] += pc;
                scores[j] += pc;
            }
        }
    }

    int best_v = -1;
    int best_score = -1;
    for (unsigned i=0; i<vertices_with_best_available_class_count.size(); i++) {
        int v = vertices_with_best_available_class_count[i];
        if (scores[i] > best_score) {
            best_score = scores[i];
            best_v = v;
        }
    }

    return best_v;
}

void expand(ColouringGraph *g, struct Solution *C,
        struct Solution *incumbent, int level, unsigned long long *expand_call_count, unsigned long long expand_call_limit,
        int num_colours, unsigned long long *available_classes_bitset,
        int *num_colours_assigned_to_vertex, int domain_num_words, int f, std::atomic_bool & terminate_early)
{
    (*expand_call_count)++;
    if (*expand_call_count >= expand_call_limit)
        return;

    if (terminate_early.load())
        return;

    if (C->size == g->n * f) {
        copy_Solution(C, incumbent);
        return;
    }

    // UNIT PROPAGATION
    int C_sz_before_unit_prop = C->size;
    std::vector<int> unit_v_stack;
    for (int i=0; i<g->n; i++) {
        int pc = bitset_popcount(available_classes_bitset + i * domain_num_words, domain_num_words);
        int num_possible_colours = pc + num_colours_assigned_to_vertex[i];
        if (pc != 0 && num_possible_colours == f) {
            unit_v_stack.push_back(i);
        } else if (num_possible_colours < f) {
            return;
        }
    }

    while (!unit_v_stack.empty()) {
        int v = unit_v_stack.back();
        unit_v_stack.pop_back();
        int colour = first_set_bit(available_classes_bitset + v * domain_num_words, domain_num_words);
        solution_colour_vtx(C, v, colour, available_classes_bitset, num_colours_assigned_to_vertex, domain_num_words, f);
        if (num_colours_assigned_to_vertex[v] != f)
            unit_v_stack.push_back(v);

        unsigned i=0;
        for (int w=0; w<g->n; w++) {
            if (i < g->adjlist[v].size() && g->adjlist[v][i] == w) {
                ++i;
                continue;
            } else if (w == v) {
                continue;
            }
            if (test_bit(available_classes_bitset + w * domain_num_words, colour)) {
                unset_bit(available_classes_bitset + w * domain_num_words, colour);
                int popcount = bitset_popcount(available_classes_bitset + w * domain_num_words, domain_num_words);
                if (popcount != 0 && popcount + num_colours_assigned_to_vertex[w] == f) {
                    unit_v_stack.push_back(w);
                } else if (popcount + num_colours_assigned_to_vertex[w] < f) {
                    solution_resize(C, C_sz_before_unit_prop);
                    return;
                }
            }
        }
    }

    if (C->size == g->n * f) {
        copy_Solution(C, incumbent);
        solution_resize(C, C_sz_before_unit_prop);
        return;
    }

    int best_v = choose_branching_vertex(g, available_classes_bitset, domain_num_words);

    std::vector<unsigned long long> colours_in_all_domains(domain_num_words, ~0ull);
    for (int i=0; i<g->n; i++)
        if (!bitset_empty(available_classes_bitset + i * domain_num_words, domain_num_words))
            bitset_intersect_with(colours_in_all_domains.data(), available_classes_bitset + i * domain_num_words, domain_num_words);

    std::vector<unsigned long long> domain_copy(domain_num_words);
    copy_bitset(available_classes_bitset + best_v * domain_num_words, domain_copy.data(), domain_num_words);

    std::vector<unsigned long long> new_available_classes_bitset(g->n * domain_num_words);
    std::vector<int> new_num_colours_assigned_to_vertex(g->n);

    bool colour_is_in_all_domains;
    do {
        int colour = first_set_bit(domain_copy.data(), domain_num_words);
        unset_bit(domain_copy.data(), colour);
        colour_is_in_all_domains = test_bit(colours_in_all_domains.data(), colour);
        
        for (int i=0; i<g->n; i++)
            new_num_colours_assigned_to_vertex[i] = num_colours_assigned_to_vertex[i];
        for (int i=0; i<g->n * domain_num_words; i++)
            new_available_classes_bitset[i] = available_classes_bitset[i];

        unsigned i=0;
        for (int w=0; w<g->n; w++) {
            if (i < g->adjlist[best_v].size() && g->adjlist[best_v][i] == w) {
                ++i;
                continue;
            } else if (w == best_v) {
                continue;
            }
            unset_bit(new_available_classes_bitset.data() + w * domain_num_words, colour);
            // We don't need to check for domain wipeout here, since any domain that was unit
            // would have been instantiated in the unit-propagation step.
        }

        solution_colour_vtx(C, best_v, colour, new_available_classes_bitset.data(), new_num_colours_assigned_to_vertex.data(), domain_num_words, f);
        expand(g, C, incumbent, level+1, expand_call_count, expand_call_limit,
                num_colours, new_available_classes_bitset.data(), new_num_colours_assigned_to_vertex.data(), domain_num_words, f, terminate_early);
        solution_pop_vtx(C);
    } while (incumbent->size < g->n * f &&
             !colour_is_in_all_domains  &&
             !bitset_empty(domain_copy.data(), domain_num_words));

    solution_resize(C, C_sz_before_unit_prop);
}

void solve(ColouringGraph *g, unsigned long long *expand_call_count, unsigned long long expand_call_limit,
        struct Solution *incumbent, int num_colours, int f, std::atomic_bool & terminate_early)
{
    struct Solution C(g->n, f);
    int domain_num_words = (num_colours + BITS_PER_WORD - 1) / BITS_PER_WORD;
    std::vector<unsigned long long> available_classes_bitset(g->n * domain_num_words);
    for (int i=0; i<g->n; i++) {
        set_first_n_bits(available_classes_bitset.data() + i * domain_num_words, num_colours);
    }
    std::vector<int> num_colours_assigned_to_vertex(g->n);
    expand(g, &C, incumbent, 0, expand_call_count, expand_call_limit,
            num_colours, available_classes_bitset.data(), num_colours_assigned_to_vertex.data(), domain_num_words, f, terminate_early);
}

// FIXME
bool is_solution_valid(ColouringGraph *g, struct Solution *solution,
        int num_colours, int f)
{
    for (int i=0; i<g->n; i++)
        for (int j=0; j<i; j++)
            if (g->adj_matrix[i][j])
                for (int k=0; k<f; k++)
                    for (int m=0; m<f; m++)
                        if (solution->vtx_colour[i][k] == solution->vtx_colour[j][m])
                            return false;
    for (int i=0; i<g->n; i++)
        for (int j=0; j<f; j++)
            if (solution->vtx_colour[i][j] < 0 || solution->vtx_colour[i][j] >= num_colours)
                return false;
    return true;
}

std::vector<int> randomised_vertex_order(const ColouringGraph & g, unsigned seed)
{
    srand(seed);

    std::vector<int> vv;
    for (int i=0; i<g.n; i++)
        vv.push_back(i);
    for (int i=g.n-1; i>=1; i--) {
        int r = rand() % (i+1);
        std::swap(vv[i], vv[r]);
    }

    return vv;
}

int find_colouring_number(const ColouringGraph & g, int f, std::atomic_bool & terminate_early)
{
    unsigned rng_seed = 0;

    std::vector<int> vv = randomised_vertex_order(g, rng_seed);
    ColouringGraph sorted_g = g.induced_subgraph(vv);

    unsigned long long expand_call_limit = 1000;
    int num_colours = 0;
    for ( ; ; num_colours++) {
        struct Solution clq(g.n, f);

        sorted_g.make_adjacency_lists();

        unsigned long long expand_call_count = 0;
        while (true) {
            if (terminate_early.load())
                return -1;

            solve(&sorted_g, &expand_call_count, expand_call_limit, &clq, num_colours, f, terminate_early);
            if (expand_call_count < expand_call_limit)
                break;
            clq.size = 0;
            expand_call_limit = expand_call_limit + expand_call_limit / 10;
            expand_call_count = 0;
            ++rng_seed;
            vv = randomised_vertex_order(g, rng_seed);
            sorted_g = g.induced_subgraph(vv);
            sorted_g.make_adjacency_lists();
        }

//        if (clq.size == sorted_g->n * arguments.fractional_level) {
//            printf("Solution");
//            struct Solution solution_for_unsorted_graph;
//            init_Solution(&solution_for_unsorted_graph, sorted_g->n);
//            solution_for_unsorted_graph.size = clq.size;
//            for (int i=0; i<sorted_g->n; i++) {
//                solution_for_unsorted_graph.vtx_colour[vv[i]] = clq.vtx_colour[i];
//            }
//            if (!is_solution_valid(g, &solution_for_unsorted_graph, num_colours))
//                fail("The solution that was found is not a colouring with the correct number of colours!!!");
//            for (int i=0; i<sorted_g->n; i++) {
//                printf(" %d", solution_for_unsorted_graph.vtx_colour[i]);
//            }
//            destroy_Solution(&solution_for_unsorted_graph);
//            printf("\n");
//        }

//        printf("%d %lld %s\n", num_colours, expand_call_count, clq.size == sorted_g.n * f ? "SATISFIABLE" : "UNSAT");

        if (clq.size == sorted_g.n * f)
            break;
    }
    return num_colours;
}

unsigned long long ColouringNumberFinder::get_search_node_count()
{
    return search_node_count;
}

int ColouringNumberFinder::get_colouring_number()
{
    return colouring_number;
}

void ColouringNumberFinder::search() {
#ifdef WITHOUT_COLOURING_UPPER_BOUND
    return;
#endif

    if (colouring_number != -1)
        return;    // solution has been found already

    std::atomic_bool terminate_early(false);

    std::vector<int> vv = randomised_vertex_order(g, rng_seed);
    ColouringGraph sorted_g = g.induced_subgraph(vv);
    sorted_g.make_adjacency_lists();

    struct Solution clq(g.n, f);

    unsigned long long local_search_node_count = 0;

    solve(&sorted_g, &local_search_node_count, local_search_node_limit, &clq, current_target_num_colours, f, terminate_early);
    search_node_count += local_search_node_count;

    if (local_search_node_count >= local_search_node_limit) {
        local_search_node_limit = local_search_node_limit + local_search_node_limit / 10;
        ++rng_seed;
        vv = randomised_vertex_order(g, rng_seed);
        sorted_g = g.induced_subgraph(vv);
        sorted_g.make_adjacency_lists();
    } else {
        if (clq.size == g.n * f) {
            colouring_number = current_target_num_colours;
        } else {
            ++current_target_num_colours;
        }
    }
}

ColouringNumberFinder::ColouringNumberFinder(const ColouringGraph & g, int f)
        : g(g), sorted_graph(0), f(f)
{
    std::vector<int> vv = randomised_vertex_order(g, rng_seed);
    ColouringGraph sorted_g = g.induced_subgraph(vv);
    sorted_g.make_adjacency_lists();
}

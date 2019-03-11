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

void add_edge(struct Graph & g, int v, int w) {
    if (!g.adj_matrix[v][w]) {
        g.adj_matrix[v][w] = true;
        g.adj_matrix[w][v] = true;
    }
}

void make_adjacency_lists(struct Graph *g)
{
    for (int i=0; i<g->n; i++)
        for (int j=0; j<g->n; j++)
            if (g->adj_matrix[i][j])
                g->adjlist[i].push_back(j);
}

struct Graph induced_subgraph(const Graph *g, std::vector<int> & vv) {
    struct Graph subg(vv.size());
    for (int i=0; i<subg.n; i++)
        for (int j=0; j<i; j++)
            if (g->adj_matrix[vv[i]][vv[j]])
                add_edge(subg, i, j);

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

void init_Solution(struct Solution *solution, int capacity)
{
    solution->vtx_colour = (int *) malloc(capacity * sizeof *solution->vtx_colour);
    solution->size = 0;
    solution->capacity = capacity;
}

void destroy_Solution(struct Solution *solution)
{
    free(solution->vtx_colour);
}

void solution_colour_vtx(struct Solution *solution, int v, int colour,
        unsigned long long *available_classes_bitset, int *num_colours_assigned_to_vertex, int domain_num_words, int f)
{
    ++solution->size;
    solution->vtx_colour[v] = colour;
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
    for (int i=0; i<src->capacity; i++)
        dest->vtx_colour[i] = src->vtx_colour[i];
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// Precondition: at least one vertex remains that can be branched on.
int choose_branching_vertex(struct Graph *original_g, unsigned long long *available_classes_bitset,
        int domain_num_words)
{
    int best_available_class_count = INT_MAX;
    std::vector<int> vertices_with_best_available_class_count;
    for (int i=0; i<original_g->n; i++) {
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
        int orig_v = vertices_with_best_available_class_count[i];
        for (unsigned j=0; j<i; j++) {
            int orig_w = vertices_with_best_available_class_count[j];
            if (!original_g->adj_matrix[orig_v][orig_w]) {
                int pc = bitset_intersection_popcount(
                        available_classes_bitset + orig_v * domain_num_words, available_classes_bitset + orig_w * domain_num_words,
                        domain_num_words);
                scores[i] += pc;
                scores[j] += pc;
            }
        }
    }

    int best_orig_v = -1;
    int best_score = -1;
    for (unsigned i=0; i<vertices_with_best_available_class_count.size(); i++) {
        int orig_v = vertices_with_best_available_class_count[i];
        if (scores[i] > best_score) {
            best_score = scores[i];
            best_orig_v = orig_v;
        }
    }

    return best_orig_v;
}

void expand(struct Graph *original_g, struct Solution *C,
        struct Solution *incumbent, int level, unsigned long long *expand_call_count, unsigned long long expand_call_limit,
        int num_colours, unsigned long long *available_classes_bitset,
        int *num_colours_assigned_to_vertex, int domain_num_words, int f)
{
    (*expand_call_count)++;
    if (*expand_call_count >= expand_call_limit)
        return;

    if (C->size == original_g->n * f) {
        copy_Solution(C, incumbent);
        return;
    }

    // UNIT PROPAGATION
    int C_sz_before_unit_prop = C->size;
    std::vector<int> unit_orig_v_stack;
    for (int i=0; i<original_g->n; i++) {
        int pc = bitset_popcount(available_classes_bitset + i * domain_num_words, domain_num_words);
        int num_possible_colours = pc + num_colours_assigned_to_vertex[i];
        if (pc != 0 && num_possible_colours == f) {
            unit_orig_v_stack.push_back(i);
        } else if (num_possible_colours < f) {
            return;
        }
    }

    while (!unit_orig_v_stack.empty()) {
        int orig_v = unit_orig_v_stack.back();
        unit_orig_v_stack.pop_back();
        int colour = first_set_bit(available_classes_bitset + orig_v * domain_num_words, domain_num_words);
        solution_colour_vtx(C, orig_v, colour, available_classes_bitset, num_colours_assigned_to_vertex, domain_num_words, f);
        if (num_colours_assigned_to_vertex[orig_v] != f)
            unit_orig_v_stack.push_back(orig_v);

        unsigned i=0;
        for (int orig_w=0; orig_w<original_g->n; orig_w++) {
            if (i < original_g->adjlist[orig_v].size() && original_g->adjlist[orig_v][i] == orig_w) {
                ++i;
                continue;
            } else if (orig_w == orig_v) {
                continue;
            }
            if (test_bit(available_classes_bitset + orig_w * domain_num_words, colour)) {
                unset_bit(available_classes_bitset + orig_w * domain_num_words, colour);
                int popcount = bitset_popcount(available_classes_bitset + orig_w * domain_num_words, domain_num_words);
                if (popcount != 0 && popcount + num_colours_assigned_to_vertex[orig_w] == f) {
                    unit_orig_v_stack.push_back(orig_w);
                } else if (popcount + num_colours_assigned_to_vertex[orig_w] < f) {
                    solution_resize(C, C_sz_before_unit_prop);
                    return;
                }
            }
        }
    }

    if (C->size == original_g->n * f) {
        copy_Solution(C, incumbent);
        solution_resize(C, C_sz_before_unit_prop);
        return;
    }

    int best_orig_v = choose_branching_vertex(original_g, available_classes_bitset, domain_num_words);

    std::vector<unsigned long long> colours_in_all_domains(domain_num_words, ~0ull);
    for (int i=0; i<original_g->n; i++)
        if (!bitset_empty(available_classes_bitset + i * domain_num_words, domain_num_words))
            bitset_intersect_with(colours_in_all_domains.data(), available_classes_bitset + i * domain_num_words, domain_num_words);

    std::vector<unsigned long long> domain_copy(domain_num_words);
    copy_bitset(available_classes_bitset + best_orig_v * domain_num_words, domain_copy.data(), domain_num_words);

    std::vector<unsigned long long> new_available_classes_bitset(original_g->n * domain_num_words);
    std::vector<int> new_num_colours_assigned_to_vertex(original_g->n);

    while (!bitset_empty(domain_copy.data(), domain_num_words)) {
        int orig_colour = first_set_bit(domain_copy.data(), domain_num_words);
        unset_bit(domain_copy.data(), orig_colour);
        bool colour_is_in_all_domains = test_bit(colours_in_all_domains.data(), orig_colour);
        
        for (int i=0; i<original_g->n; i++)
            new_num_colours_assigned_to_vertex[i] = num_colours_assigned_to_vertex[i];
        for (int i=0; i<original_g->n * domain_num_words; i++)
            new_available_classes_bitset[i] = available_classes_bitset[i];

        unsigned i=0;
        for (int orig_w=0; orig_w<original_g->n; orig_w++) {
            if (i < original_g->adjlist[best_orig_v].size() && original_g->adjlist[best_orig_v][i] == orig_w) {
                ++i;
                continue;
            } else if (orig_w == best_orig_v) {
                continue;
            }
            unset_bit(new_available_classes_bitset.data() + orig_w * domain_num_words, orig_colour);
            // We don't need to check for domain wipeout here, since any domain that was unit
            // would have been instantiated in the unit-propagation step.
        }

        solution_colour_vtx(C, best_orig_v, orig_colour, new_available_classes_bitset.data(), new_num_colours_assigned_to_vertex.data(), domain_num_words, f);
        expand(original_g, C, incumbent, level+1, expand_call_count, expand_call_limit,
                num_colours, new_available_classes_bitset.data(), new_num_colours_assigned_to_vertex.data(), domain_num_words, f);
        solution_pop_vtx(C);

        if (incumbent->size == original_g->n * f)
            break;

        if (colour_is_in_all_domains)
            break;
    }

    solution_resize(C, C_sz_before_unit_prop);
}

void solve(struct Graph *original_g, unsigned long long *expand_call_count, unsigned long long expand_call_limit,
        struct Solution *incumbent, int num_colours, int f)
{
    struct Solution C;
    init_Solution(&C, original_g->n);
    int domain_num_words = (num_colours + BITS_PER_WORD - 1) / BITS_PER_WORD;
    std::vector<unsigned long long> available_classes_bitset(original_g->n * domain_num_words);
    for (int i=0; i<original_g->n; i++) {
        set_first_n_bits(available_classes_bitset.data() + i * domain_num_words, num_colours);
    }
    std::vector<int> num_colours_assigned_to_vertex(original_g->n);
    expand(original_g, &C, incumbent, 0, expand_call_count, expand_call_limit,
            num_colours, available_classes_bitset.data(), num_colours_assigned_to_vertex.data(), domain_num_words, f);
    destroy_Solution(&C);
}

// FIXME
bool is_solution_valid(struct Graph *original_g, struct Solution *solution,
        int num_colours)
{
    for (int i=0; i<original_g->n; i++)
        for (int j=0; j<i; j++)
            if (original_g->adj_matrix[i][j] && solution->vtx_colour[i] == solution->vtx_colour[j])
                return false;
    for (int i=0; i<original_g->n; i++)
        if (solution->vtx_colour[i] < 0 || solution->vtx_colour[i] >= num_colours)
            return false;
    return true;
}

std::vector<int> randomised_vertex_order(const Graph & g, unsigned seed)
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

int find_colouring_number(const Graph & g, int f)
{
    unsigned rng_seed = 0;

    std::vector<int> vv = randomised_vertex_order(g, rng_seed);
    struct Graph sorted_g = induced_subgraph(&g, vv);

    unsigned long long expand_call_limit = 1000;
    int num_colours = 0;
    for ( ; ; num_colours++) {
        struct Solution clq;
        init_Solution(&clq, sorted_g.n);

        make_adjacency_lists(&sorted_g);

        unsigned long long expand_call_count = 0;
        while (true) {
            solve(&sorted_g, &expand_call_count, expand_call_limit, &clq, num_colours, f);
            if (expand_call_count < expand_call_limit)
                break;
            clq.size = 0;
            expand_call_limit = expand_call_limit + expand_call_limit / 10;
            expand_call_count = 0;
            ++rng_seed;
            vv = randomised_vertex_order(g, rng_seed);
            sorted_g = induced_subgraph(&g, vv);
            make_adjacency_lists(&sorted_g);
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

        printf("%d %lld %s\n", num_colours, expand_call_count, clq.size == sorted_g.n * f ? "SATISFIABLE" : "UNSAT");

        destroy_Solution(&clq);

        if (clq.size == sorted_g.n * f)
            break;
    }
    return num_colours;
}


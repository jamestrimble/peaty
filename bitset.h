#ifndef BITSET_H
#define BITSET_H

#include "graph.h"

static void set_first_n_bits(vector<unsigned long long> & bitset, int n)
{
    int i = 0;
    while (n >= 64) {
        bitset[i] = ~0ull;
        n -= 64;
        ++i;
    }
    if (n > 0)
        bitset[i] = (1ull << n) - 1;
}

static bool test_bit(const vector<unsigned long long> & bitset, int bit)
{
    return 0 != (bitset[bit/BITS_PER_WORD] & (1ull << (bit%BITS_PER_WORD)));
}

static void set_bit(vector<unsigned long long> & bitset, int bit)
{
    bitset[bit/BITS_PER_WORD] |= (1ull << (bit%BITS_PER_WORD));
}

static void unset_bit(vector<unsigned long long> & bitset, int bit)
{
    bitset[bit/BITS_PER_WORD] &= ~(1ull << (bit%BITS_PER_WORD));
}

static void unset_bit_if(vector<unsigned long long> & bitset, int bit, bool condition)
{
    bitset[bit/BITS_PER_WORD] &= ~((unsigned long long) condition << (bit%BITS_PER_WORD));
}

static int bitset_popcount(const vector<unsigned long long> & bitset, int num_words)
{
    int count = 0;
    for (int i=0; i<num_words; i++)
        count += __builtin_popcountll(bitset[i]);
    return count;
}

static int bitset_intersection_popcount(const vector<unsigned long long> & bitset1,
        const vector<unsigned long long> & bitset2, int num_words)
{
    int count = 0;
    for (int i=0; i<num_words; i++)
        count += __builtin_popcountll(bitset1[i] & bitset2[i]);
    return count;
}

static bool bitset_empty(const vector<unsigned long long> & bitset, int num_words)
{
    for (int i=0; i<num_words; i++)
        if (bitset[i] != 0)
            return false;
    return true;
}

static int first_set_bit(const vector<unsigned long long> & bitset,
                         int num_words)
{
    for (int i=0; i<num_words; i++)
        if (bitset[i] != 0)
            return i*BITS_PER_WORD + __builtin_ctzll(bitset[i]);
    return -1;
}

static int last_set_bit(const vector<unsigned long long> & bitset,
                         int num_words)
{
    for (int i=num_words; i--; )
        if (bitset[i] != 0)
            return i*BITS_PER_WORD + BITS_PER_WORD - 1 - __builtin_clzll(bitset[i]);
    return -1;
}

static bool have_non_empty_intersection(const vector<unsigned long long> & bitset1,
                                     const vector<unsigned long long> & bitset2,
                                     int num_words)
{
    for (int i=0; i<num_words; i++)
        if (bitset1[i] & bitset2[i])
            return true;
    return false;
}

static bool have_empty_intersection(const vector<unsigned long long> & bitset1,
                                     const vector<unsigned long long> & bitset2,
                                     int num_words)
{
    for (int i=0; i<num_words; i++)
        if (bitset1[i] & bitset2[i])
            return false;
    return true;
}

static int first_nonzero_in_intersection(const vector<unsigned long long> & bitset1,
                                     const vector<unsigned long long> & bitset2,
                                     int num_words)
{
    for (int i=0; i<num_words; i++) {
        unsigned long long word_intersection = bitset1[i] & bitset2[i];
        if (word_intersection)
            return i*BITS_PER_WORD + __builtin_ctzll(word_intersection);
    }
    return -1;
}

static void bitset_intersection(const vector<unsigned long long> & src1,
                                     const vector<unsigned long long> & src2,
                                     vector<unsigned long long> & dst,
                                     int num_words)
{
    for (int i=0; i<num_words; i++)
        dst[i] = src1[i] & src2[i];
}

static void bitset_intersect_with(vector<unsigned long long> & bitset1,
                                     const vector<unsigned long long> & bitset2,
                                     int num_words)
{
    for (int i=0; i<num_words; i++)
        bitset1[i] &= bitset2[i];
}

static void bitset_intersect_with_complement(vector<unsigned long long> & bitset,
                                     const vector<unsigned long long> & bitset2,
                                     int num_words)
{
    for (int i=0; i<num_words; i++)
        bitset[i] &= ~bitset2[i];
}

static void bitset_intersection_with_complement(const vector<unsigned long long> & src1,
                                     const vector<unsigned long long> & src2,
                                     vector<unsigned long long> & dst,
                                     int num_words)
{
    for (int i=0; i<num_words; i++)
        dst[i] = src1[i] & ~src2[i];
}

static void copy_bitset(const vector<unsigned long long> & src,
                        vector<unsigned long long> & dest,
                        int num_words)
{
    for (int i=0; i<num_words; i++)
        dest[i] = src[i];
}

static void clear_bitset(vector<unsigned long long> & bitset,
                        int num_words)
{
    for (int i=0; i<num_words; i++)
        bitset[i] = 0ull;
}

static int calc_numwords(const vector<unsigned long long> & bitset, int graph_num_words)
{
    for (int i=graph_num_words; i--; )
        if (bitset[i] != 0)
            return i+1;
    return 0;
}

template<typename F>
static void bitset_foreach(const vector<unsigned long long> & bitset, F f, int numwords)
{
        for (int i=0; i<numwords; i++) {
            unsigned long long word = bitset[i];
            while (word) {
                int bit = __builtin_ctzll(word);
                word ^= (1ull << bit);
                int v = i*BITS_PER_WORD + bit;
                f(v);
            }
        }
}

#endif

#ifndef FAST_INT_QUEUE_H
#define FAST_INT_QUEUE_H

#include <assert.h>

// This implementation is very restrictive: it assumes that there are
// at most `capacity` calls to enqueue() after each call to clear()
struct FastIntQueue {
    vector<int> vals;
    int start;
    int end;
    FastIntQueue(int capacity)
            : vals(capacity), start(0), end(0)
    {
    }
    auto enqueue(int val) -> void
    {
        vals[end++] = val;
    }
    auto dequeue() -> int
    {
        assert (start != end);
        return vals[start++];
    }
    auto clear() -> void
    {
        start = 0;
        end = 0;
    }
    auto empty() -> bool
    {
        return start == end;
    }
};

#endif

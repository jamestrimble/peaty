#ifndef INT_STACK_WITHOUT_DUPS_H
#define INT_STACK_WITHOUT_DUPS_H

struct IntStackWithoutDups {
    vector<int> vals;
    vector<unsigned char> on_stack;  // really a vector of boolean values
    IntStackWithoutDups(int max_size)
            : on_stack(max_size)
    {
        vals.reserve(max_size);
    }
    auto push(int val) -> void
    {
        if (!on_stack[val]) {
            vals.push_back(val);
            on_stack[val] = 1;
        }
    }
    auto clear() -> void
    {
        for (int val : vals)
            on_stack[val] = 0;
        vals.clear();
    }
    auto empty() -> bool
    {
        return vals.empty();
    }
};

#endif

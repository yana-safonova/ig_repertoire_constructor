#pragma once
#include <cassert>

namespace core {
    template<typename Titer, typename Tf>
    static auto max_map(Titer b, Titer e, const Tf &f) -> decltype(f(*b)) { // TODO Add decay
        assert(b != e);

        auto result = f(*b++);
        for(;b != e; ++b) {
            auto cur = f(*b);
            if (result <  cur) {
                result = cur;
            }
        }

        return result;
    }
}
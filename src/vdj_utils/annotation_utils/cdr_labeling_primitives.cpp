#include <verify.hpp>

#include "cdr_labeling_primitives.hpp"

namespace annotation_utils {
    size_t CDRRange::length() const  {
        VERIFY_MSG(Full(), "Start pos (" << start_pos << ") or end pos (" << end_pos << ") is not defined");
        return end_pos - start_pos + 1;
    }

    std::ostream& operator<<(std::ostream& out, const CDRRange &obj) {
        out << "Start pos: " << obj.start_pos << ", end pos: " << obj.end_pos;
        return out;
    }
}
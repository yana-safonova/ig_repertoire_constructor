#include <verify.hpp>

#include "cdr_labeling_primitives.hpp"

namespace annotation_utils {
    size_t CDRRange::length() const  {
        VERIFY_MSG(Full(), "Start pos (" << start_pos << ") or end pos (" << end_pos << ") is not defined");
        return end_pos - start_pos + 1;
    }
}
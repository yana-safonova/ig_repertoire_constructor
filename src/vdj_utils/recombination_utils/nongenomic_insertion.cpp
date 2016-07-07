#include "nongenomic_insertion.hpp"

namespace recombination_utils {

std::ostream &operator<<(std::ostream &out, const NongenomicInsertion &insertion) {
    out << "Start: " << insertion.StartPosition() << ", end: " << insertion.EndPosition() << ", len: " <<
        insertion.length() << ", valid: " << insertion.Valid();
    return out;
}

}
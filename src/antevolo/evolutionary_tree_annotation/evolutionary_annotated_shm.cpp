#include "evolutionary_annotated_shm.hpp"

namespace antevolo {
    std::ostream& operator<<(std::ostream &out, const EvolutionaryAnnotatedSHM &shm) {
        out << shm.shm << ", synonymous: " << shm.synonymous << ", fixed: " << shm.fixed;
        return out;
    }
}
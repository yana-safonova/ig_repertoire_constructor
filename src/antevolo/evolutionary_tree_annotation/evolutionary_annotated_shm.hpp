#pragma once

#include <annotation_utils/shm_annotation/shm_annotation.hpp>

namespace antevolo {
    struct EvolutionaryAnnotatedSHM {
        annotation_utils::SHM shm; // base SHM
        bool synonymous; // true if shm does not change amino acid in parent node
        bool fixed; // fixed if it is presented in more than one node in a clonal tree

        EvolutionaryAnnotatedSHM(annotation_utils::SHM shm, bool synonymous, bool fixed) :
                shm(shm), synonymous(synonymous), fixed(fixed) { }
    };

    std::ostream& operator<<(std::ostream &out, const EvolutionaryAnnotatedSHM &shm);
}
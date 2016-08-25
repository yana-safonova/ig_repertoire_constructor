#pragma once

#include <annotation_utils/shm_annotation/shm_annotation.hpp>

namespace antevolo {
    struct EvolutionaryAnnotatedSHM {
        annotation_utils::SHM shm; // base SHM
        bool synonymous; // true if shm does not change amino acid in parent node
        bool fixed; // fixed if it is presented in more than one node in a clonal tree

        EvolutionaryAnnotatedSHM(annotation_utils::SHM shm, bool synonymous, bool fixed) :
                shm(shm), synonymous(synonymous), fixed(fixed) { }

        bool operator==(const EvolutionaryAnnotatedSHM& annotated_shm) const {
            return shm == annotated_shm.shm;
        }
    };

    struct EvolutionaryAnnotatedSHMComparator {
        bool operator()(const EvolutionaryAnnotatedSHM &shm1, const EvolutionaryAnnotatedSHM& shm2) {
            return annotation_utils::TrivialSHMComparator()(shm1.shm, shm2.shm);
        }
    };

    std::ostream& operator<<(std::ostream &out, const EvolutionaryAnnotatedSHM &shm);

    //----------------------------------------------------------------------------------
    struct CDR3SHM {
        size_t src_pos;
        size_t dst_pos;
        size_t rel_src_pos;
        size_t rel_dst_pos;
        char src_nucl;
        char dst_nucl;
        char src_aa;
        char dst_aa;

        CDR3SHM(size_t src_pos, size_t dst_pos,
                size_t rel_src_pos, size_t rel_dst_pos,
                char src_nucl, char dst_nucl,
                char src_aa, char dst_aa) :
                src_pos(src_pos), dst_pos(dst_pos),
                rel_src_pos(rel_src_pos), rel_dst_pos(rel_dst_pos),
                src_nucl(src_nucl), dst_nucl(dst_nucl),
                src_aa(src_aa), dst_aa(dst_aa) { }

        bool IsSynonymous() const { return src_aa == dst_aa; }

        bool operator==(const CDR3SHM& cdr3_shm) {
            return rel_src_pos == cdr3_shm.rel_src_pos and rel_dst_pos == cdr3_shm.rel_dst_pos and
                    src_nucl == cdr3_shm.src_nucl and dst_nucl == cdr3_shm.dst_nucl;
        }
    };

    struct CDR3SHMComparator {
        bool operator()(const CDR3SHM& shm1, const CDR3SHM& shm2);
    };
}
#pragma once

#include <verify.hpp>
#include "../germline_utils/germline_gene_type.hpp"

namespace annotation_utils {
    enum SHMType { UnknownSHM, SubstitutionSHM, InsertionSHM, DeletionSHM };

    struct SHM {
        SHMType shm_type;
        size_t gene_pos;
        size_t read_pos;
        char gene_nucl;
        char read_nucl;
        char gene_aa;
        char read_aa;

        void ComputeType();

        SHM(size_t gene_pos, size_t read_pos, char gene_nucl,
            char read_nucl, char gene_aa, char read_aa) : gene_pos(gene_pos), read_pos(read_pos),
                                                          gene_nucl(gene_nucl), read_nucl(read_nucl), gene_aa(gene_aa),
                                                          read_aa(read_aa) {
            ComputeType();
        }

        bool IsSynonymous() const {
            VERIFY_MSG(shm_type == SHMType::SubstitutionSHM, "SHM is not substitution");
            return gene_aa == read_aa;
        }

        bool ToStopCodon() const {
            VERIFY_MSG(shm_type == SHMType::SubstitutionSHM, "SHM is not substitution");
            return read_aa == '*';
        }
    };

    std::ostream& operator<<(std::ostream &out, const SHM& shm);

    // class stores SHMs in the order of increasing positions
    class GeneSegmentSHMs {
        germline_utils::SegmentType segment_type_;
        std::vector<SHM> shms_;

        void CheckConsistencyFatal(SHM shm);

    public:
        GeneSegmentSHMs(germline_utils::SegmentType segment_type) : segment_type_(segment_type) { }

        void AddSHM(SHM shm);

        size_t size() const { return shms_.size(); }

        const SHM& operator[](size_t index) const;

        typedef std::vector<SHM>::const_iterator SHMConstIterator;

        SHMConstIterator cbegin() const { return shms_.cbegin(); }

        SHMConstIterator cend() const { return shms_.cend(); }

        bool operator==(const SHM& obj) const;

        germline_utils::SegmentType SegmentType() const { return segment_type_; }
    };

    std::ostream& operator<<(std::ostream &out, const GeneSegmentSHMs& shms);
}
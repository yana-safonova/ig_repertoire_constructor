#pragma once

#include <verify.hpp>
#include <read_archive.hpp>
#include <germline_utils/germline_databases/immune_gene_database.hpp>

namespace annotation_utils {
    enum SHMType { UnknownSHM, SubstitutionSHM, InsertionSHM, DeletionSHM };

    std::ostream& operator<<(std::ostream& out, const SHMType &shm_type);

    struct SHM {
        germline_utils::SegmentType segment_type;
        SHMType shm_type;
        size_t gene_nucl_pos;
        size_t read_nucl_pos;
        char gene_nucl;
        char read_nucl;
        char gene_aa;
        char read_aa;

        void ComputeType();

        SHM(germline_utils::SegmentType segment_type,
            size_t gene_nucl_pos,
            size_t read_nucl_pos,
            char gene_nucl,
            char read_nucl,
            char gene_aa,
            char read_aa) : segment_type(segment_type),
                            gene_nucl_pos(gene_nucl_pos),
                            read_nucl_pos(read_nucl_pos),
                            gene_nucl(gene_nucl),
                            read_nucl(read_nucl),
                            gene_aa(gene_aa),
                            read_aa(read_aa) {
            ComputeType();
        }

        bool IsSynonymous() const {
            if(shm_type != SHMType::SubstitutionSHM)
                return false;
            return gene_aa == read_aa;
        }

        bool ToStopCodon() const {
            if(shm_type != SHMType::SubstitutionSHM)
                return false;
            return read_aa == '*';
        }

        bool operator==(const SHM& shm) const;

        bool operator!=(const SHM& shm) const;

        void AppendInMixcrFormat(std::ostream& out) const;
    };

    bool operator<(const SHM &left, const SHM &right);

    struct TrivialSHMComparator {
        bool operator()(const SHM& shm1, const SHM& shm2);
    };

    std::ostream& operator<<(std::ostream &out, const SHM& shm);

    // class stores SHMs in the order of increasing positions
    class GeneSegmentSHMs {
        core::Read read_;
        const germline_utils::ImmuneGene *immune_gene_;

        std::vector<SHM> shms_;

        void CheckConsistencyFatal(SHM shm);

    public:
        GeneSegmentSHMs(core::Read read,
                        const germline_utils::ImmuneGene &immune_gene) :
                read_(read),
                immune_gene_(&immune_gene){ }

        void AddSHM(SHM shm);

        size_t size() const { return shms_.size(); }

        const SHM& operator[](size_t index) const;

        typedef std::vector<SHM>::const_iterator SHMConstIterator;

        SHMConstIterator cbegin() const { return shms_.cbegin(); }

        SHMConstIterator cend() const { return shms_.cend(); }

        bool operator==(const SHM& obj) const;

        germline_utils::SegmentType SegmentType() const { return immune_gene_->Segment(); }

        const core::Read& Read() const { return read_; }

        const germline_utils::ImmuneGene& ImmuneGene() const { return *immune_gene_; }

        /**
         * Appends MiXCR-like alignment.
         * See <a href="http://mixcr.readthedocs.io/en/latest/appendix.html#ref-encoding">MiXCR documentation</a> for description.
         */
        void AppendInMixcrFormat(std::ostream& out) const;
    };

    std::ostream& operator<<(std::ostream &out, const GeneSegmentSHMs& shms);
}
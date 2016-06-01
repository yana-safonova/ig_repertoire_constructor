#pragma once

#include <germline_utils/germline_databases/custom_gene_database.hpp>

namespace cdr_labeler {
    // start and end positions are inclusive
    struct CDRRange {
        size_t start_pos;
        size_t end_pos;

        CDRRange() :
                start_pos(size_t(-1)),
                end_pos(size_t(-1)) { }

        CDRRange(size_t start_pos, size_t end_pos) :
                start_pos(start_pos),
                end_pos(end_pos) {
            if(start_pos != size_t(-1))
                VERIFY_MSG(start_pos < end_pos, "Start position (" << start_pos <<
                    ") exceeds end position (" << end_pos << ")");
        }

        bool Valid() const {
            if(start_pos == size_t(-1) and end_pos == size_t(-1))
                return false;
            if(start_pos != size_t(-1) and end_pos != size_t(-1))
                return start_pos < end_pos;
            return true;
        }

        bool Empty() const { return !Valid(); }
    };

    struct CDRLabeling {
        CDRRange cdr1;
        CDRRange cdr2;
        CDRRange cdr3;

        CDRLabeling() : cdr1(), cdr2(), cdr3() { }

        CDRLabeling(CDRRange cdr1, CDRRange cdr2, CDRRange cdr3) :
                cdr1(cdr1), cdr2(cdr2), cdr3(cdr3) { }

        bool Empty() const { return cdr1.Empty() and cdr2.Empty() and cdr3.Empty(); }
    };

    class DbCDRLabeling {
        std::vector<CDRLabeling> cdr_labelings_;
        std::unordered_map<std::string, size_t> gene_name_index_map_;

        size_t num_empty_labelings_;

    public:
        DbCDRLabeling() : num_empty_labelings_() { }

        void AddGeneLabeling(const germline_utils::ImmuneGene &immune_gene, CDRLabeling labeling);

        CDRLabeling GetLabelingByGene(const germline_utils::ImmuneGene &immune_gene) const;

        bool CDRLabelingIsEmpty(const germline_utils::ImmuneGene &immune_gene) const;

        size_t NumEmptyLabelings() const { return num_empty_labelings_; }
    };
}
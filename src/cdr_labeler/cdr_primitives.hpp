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
            return start_pos < end_pos;
        }
    };

    struct CDRLabeling {
        CDRRange cdr1;
        CDRRange cdr2;
        CDRRange cdr3;

        CDRLabeling(CDRRange cdr1, CDRRange cdr2, CDRRange cdr3) :
                cdr1(cdr1), cdr2(cdr2), cdr3(cdr3) { }
    };

    // todo: make it template
    class DbCDRLabeling {
        const germline_utils::CustomGeneDatabase &gene_db_;
        std::vector<CDRLabeling> cdr_labelings_;

    public:
        DbCDRLabeling(const germline_utils::CustomGeneDatabase &gene_db) : gene_db_(gene_db) { }
    };
}
//
// Created by Andrew Bzikadze on 3/16/17.
//

#pragma once

#include <tuple>
#include <memory>

#include "germline_utils/germline_databases/custom_gene_database.hpp"

namespace ig_simulator {

using VDJ_GenesIndexTuple = std::tuple<size_t, size_t, size_t>;

class AbstractVDJGeneChooser {
protected:
    const germline_utils::CustomGeneDatabase *v_db_p_;
    const germline_utils::CustomGeneDatabase *d_db_p_;
    const germline_utils::CustomGeneDatabase *j_db_p_;

    // This variable defines whether D segment is generated
    bool is_vdj;

private:
    AbstractVDJGeneChooser(const germline_utils::CustomGeneDatabase *v_db_p,
                           const germline_utils::CustomGeneDatabase *d_db_p,
                           const germline_utils::CustomGeneDatabase *j_db_p,
                           bool is_vdj) :
        v_db_p_(v_db_p), d_db_p_(d_db_p), j_db_p_(j_db_p), is_vdj(is_vdj)
    { }

public:
    AbstractVDJGeneChooser(const germline_utils::CustomGeneDatabase &v_db,
                           const germline_utils::CustomGeneDatabase &d_db,
                           const germline_utils::CustomGeneDatabase &j_db) :
        AbstractVDJGeneChooser(&v_db, &d_db, &j_db, true)
    {
        VERIFY(v_db.size() > 0);
        VERIFY(d_db.size() > 0);
        VERIFY(j_db.size() > 0);
    }

    AbstractVDJGeneChooser(const germline_utils::CustomGeneDatabase &v_db,
                           const germline_utils::CustomGeneDatabase &j_db) :
        AbstractVDJGeneChooser(&v_db, nullptr, &j_db, false)
    {
        VERIFY(v_db.size() > 0);
        VERIFY(j_db.size() > 0);
    }

    virtual VDJ_GenesIndexTuple ChooseGenes() const = 0;

    /**
     * This method suggests whether D segment is generated.
     * If `false` then second component of `VDJ_GenesIndexTuple`
     * returned by `ChooseGenes()` will be size_t(-1).
     */
    bool IsVDJ() const { return is_vdj; }

    virtual ~AbstractVDJGeneChooser() { };
};

using AbstractVDJGeneChooserPtr = std::unique_ptr<AbstractVDJGeneChooser>;
} // End namespace ig_simulator

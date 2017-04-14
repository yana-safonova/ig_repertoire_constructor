//
// Created by Andrew Bzikadze on 3/16/17.
//

#pragma once

#include <tuple>
#include <memory>

#include "germline_utils/germline_databases/custom_gene_database.hpp"
#include "ig_simulator_utils.hpp"

namespace ig_simulator {

using VDJ_GenesIndexTuple = std::tuple<size_t, size_t, size_t>;

class AbstractVDJGeneChooser {
protected:
    const germline_utils::CustomGeneDatabase *v_db_p_;
    const germline_utils::CustomGeneDatabase *d_db_p_;
    const germline_utils::CustomGeneDatabase *j_db_p_;

    // This variable defines whether D segment is generated
    // d_dp_p_ MUST be nullptr if is_vdj == false
    bool is_vdj;

public:
    AbstractVDJGeneChooser(const std::vector<germline_utils::CustomGeneDatabase>& db):
            v_db_p_(&db.front()),
            d_db_p_(nullptr),
            j_db_p_(&db.back()),
            is_vdj(false)
    {
        VERIFY(db.size() >= 2 and db.size() <= 3);

        if (db.size() == 3) {
            d_db_p_ = &db[1];
            is_vdj = true;
        }
    }

    virtual VDJ_GenesIndexTuple ChooseGenes() const = 0;

    /**
     * This method suggests whether D segment is generated.
     * If `false` then second component of `VDJ_GenesIndexTuple`
     * returned by `ChooseGenes()` will be size_t(-1).
     */
    bool IsVDJ() const {
        if (not is_vdj) {
            VERIFY(d_db_p_ != nullptr);
        }
        return is_vdj;
    }

    virtual ~AbstractVDJGeneChooser() { };
};

using AbstractVDJGeneChooserCPtr = std::unique_ptr<const AbstractVDJGeneChooser>;

} // End namespace ig_simulator

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

public:
    AbstractVDJGeneChooser(const germline_utils::CustomGeneDatabase *v_db_p = nullptr,
                           const germline_utils::CustomGeneDatabase *d_db_p = nullptr,
                           const germline_utils::CustomGeneDatabase *j_db_p = nullptr) :
        v_db_p_(v_db_p), d_db_p_(d_db_p), j_db_p_(j_db_p)
    { }

    AbstractVDJGeneChooser(const germline_utils::CustomGeneDatabase &v_db,
                           const germline_utils::CustomGeneDatabase &d_db,
                           const germline_utils::CustomGeneDatabase &j_db) :
            v_db_p_(&v_db), d_db_p_(&d_db), j_db_p_(&j_db) {
        if (d_db.size() == 0) {
            d_db_p_ = nullptr;
        }
    }

    virtual VDJ_GenesIndexTuple ChooseGenes() const = 0;

    virtual ~AbstractVDJGeneChooser() { };
};

using AbstractVDJGeneChooserPtr = std::unique_ptr<AbstractVDJGeneChooser>;
} // End namespace ig_simulator

//
// Created by Andrew Bzikadze on 3/16/17.
//

#pragma once

#include <tuple>

#include "germline_utils/germline_databases/custom_gene_database.hpp"

class AbstractVDJGeneChooser {
private:
    const germline_utils::CustomGeneDatabase *v_db_p_;
    const germline_utils::CustomGeneDatabase *d_db_p_;
    const germline_utils::CustomGeneDatabase *j_db_p_;

public:
    using VDJ_GenesIndexTuple = std::tuple<size_t, size_t, size_t>;

public:
    AbstractVDJGeneChooser(const germline_utils::CustomGeneDatabase *v_db_p = nullptr,
                           const germline_utils::CustomGeneDatabase *d_db_p = nullptr,
                           const germline_utils::CustomGeneDatabase *j_db_p = nullptr) :
        v_db_p_(v_db_p), d_db_p_(d_db_p), j_db_p_(j_db_p)
    { }

    virtual VDJ_GenesIndexTuple ChooseGenes() = 0;
};


class AbstractVJGeneChooser {
private:
    const germline_utils::CustomGeneDatabase *v_db_p_;
    const germline_utils::CustomGeneDatabase *j_db_p_;

public:
    using VJ_GenesIndexTuple = std::tuple<size_t, size_t>;

public:
    AbstractVJGeneChooser(const germline_utils::CustomGeneDatabase *v_db_p = nullptr,
                          const germline_utils::CustomGeneDatabase *j_db_p = nullptr) :
        v_db_p_(v_db_p), j_db_p_(j_db_p)
    { }

    virtual VJ_GenesIndexTuple ChooseGenes() = 0;
};

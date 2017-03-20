//
// Created by Andrew Bzikadze on 3/17/17.
//

#pragma once

#include <random>

#include "abstract_gene_chooser.hpp"

namespace ig_simulator {

class UniformVDJGeneChooser : public AbstractVDJGeneChooser {
public:
    UniformVDJGeneChooser(const germline_utils::CustomGeneDatabase *v_db_p = nullptr,
                          const germline_utils::CustomGeneDatabase *d_db_p = nullptr,
                          const germline_utils::CustomGeneDatabase *j_db_p = nullptr) :
            AbstractVDJGeneChooser(v_db_p, d_db_p, j_db_p)
    { }

    UniformVDJGeneChooser(const germline_utils::CustomGeneDatabase &v_db,
                          const germline_utils::CustomGeneDatabase &d_db,
                          const germline_utils::CustomGeneDatabase &j_db) :
        AbstractVDJGeneChooser(v_db, d_db, j_db)
    { }

    virtual VDJ_GenesIndexTuple ChooseGenes() const override;

    virtual ~UniformVDJGeneChooser() { }
};

} // End namespace ig_simulator

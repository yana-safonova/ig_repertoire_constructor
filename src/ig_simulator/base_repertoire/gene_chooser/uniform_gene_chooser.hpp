//
// Created by Andrew Bzikadze on 3/17/17.
//

#pragma once

#include <random>

#include "abstract_gene_chooser.hpp"

namespace ig_simulator {

class UniformVDJGeneChooser : public AbstractVDJGeneChooser {
public:
    explicit UniformVDJGeneChooser(const std::vector<germline_utils::CustomGeneDatabase>& db):
        AbstractVDJGeneChooser(db)
    { }

    VDJ_GenesIndexTuple ChooseGenes() const override;
};

} // End namespace ig_simulator

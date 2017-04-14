//
// Created by Andrew Bzikadze on 3/31/17.
//

#pragma once

#include "abstract_gene_chooser.hpp"
#include "ig_simulator_config.hpp"

namespace ig_simulator {

AbstractVDJGeneChooserCPtr
get_gene_chooser(const GeneChooserParams& config,
                 const std::vector<germline_utils::CustomGeneDatabase>& db);

} // End namespace ig_simulator

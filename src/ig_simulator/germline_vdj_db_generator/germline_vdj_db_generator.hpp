//
// Created by Andrew Bzikadze on 3/15/17.
//

#pragma once

#include "germline_db_generator.hpp"

namespace ig_simulator {

class GermlineVDJDbGenerator {
    const VJFinderConfig::IOParams::InputParams::GermlineInput &germ_input_;
    const VJFinderConfig::AlgorithmParams::GermlineParams &germ_params_;

    std::vector<germline_utils::ChainType> chain_types_;
    std::vector<std::string> v_genes_fnames_;
    std::vector<std::string> j_genes_fnames_;

    void GenerateGeneFnames();

public:
    GermlineDbGenerator(const VJFinderConfig::IOParams::InputParams::GermlineInput &germ_input,
                        const VJFinderConfig::AlgorithmParams::GermlineParams &germ_params) :
        germ_input_(germ_input),
        germ_params_(germ_params) {
        GenerateGeneFnames();
    }

    germline_utils::CustomGeneDatabase GenerateVariableDb();

    germline_utils::CustomGeneDatabase GenerateJoinDb();
};


} // End namespace ig_simulator
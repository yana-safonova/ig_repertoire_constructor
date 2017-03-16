#pragma once

#include "germline_utils/germline_config.hpp"
#include <germline_utils/germline_databases/custom_gene_database.hpp>

namespace vj_finder {
    class GermlineDbGenerator {
        const germline_utils::GermlineInput &germ_input_;
        const germline_utils::GermlineParams &germ_params_;

        std::vector<germline_utils::ChainType> chain_types_;
        std::vector<std::string> v_genes_fnames_;
        std::vector<std::string> j_genes_fnames_;

        void GenerateGeneFnames();

    public:
        GermlineDbGenerator(const germline_utils::GermlineInput &germ_input,
                            const germline_utils::GermlineParams &germ_params) :
                germ_input_(germ_input),
                germ_params_(germ_params) {
            GenerateGeneFnames();
        }

        germline_utils::CustomGeneDatabase GenerateVariableDb();

        germline_utils::CustomGeneDatabase GenerateJoinDb();
    };
}
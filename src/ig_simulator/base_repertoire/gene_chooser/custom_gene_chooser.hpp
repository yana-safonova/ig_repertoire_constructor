//
// Created by Andrew Bzikadze on 4/26/17.
//

#pragma once

#include <random>
#include "random_generator.hpp"
#include "abstract_gene_chooser.hpp"
#include <boost/tokenizer.hpp>
#include <ig_simulator_config.hpp>

namespace ig_simulator {

class CustomGeneChooser final : public AbstractVDJGeneChooser {
private:
    mutable std::discrete_distribution<size_t> v_distr;
    mutable std::discrete_distribution<size_t> d_distr;
    mutable std::discrete_distribution<size_t> j_distr;

private:
    static std::vector<double> ReadProbabilities(const std::string& filename,
                                                 const germline_utils::CustomGeneDatabase& db);

    static std::discrete_distribution<size_t> GetDistr(const std::string& filename,
                                                       const germline_utils::CustomGeneDatabase& db);

public:
    CustomGeneChooser(const std::vector<germline_utils::CustomGeneDatabase>& db,
                      const std::string& v_genes_probs,
                      const std::string& d_genes_probs,
                      const std::string& j_genes_probs);

    CustomGeneChooser(const std::vector<germline_utils::CustomGeneDatabase>& db,
                      const GeneChooserParams::CustomGeneChooserParams& config);

    VDJ_GenesIndexTuple ChooseGenes() const override;
};

} // End namespace ig_simulator

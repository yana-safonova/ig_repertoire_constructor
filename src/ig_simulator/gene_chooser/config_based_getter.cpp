//
// Created by Andrew Bzikadze on 3/31/17.
//

#include "config_based_getter.hpp"
#include "uniform_gene_chooser.hpp"


namespace ig_simulator {

AbstractVDJGeneChooserPtr get_gene_chooser(const GeneChooserParams& config,
                                           const std::vector<const germline_utils::CustomGeneDatabase *>& db)
{
    if (config.method == GeneChooserMethod::Uniform)
        return AbstractVDJGeneChooserPtr(new UniformVDJGeneChooser(db));
    VERIFY(false);
}

} // End namespace ig_simulator


//
// Created by Andrew Bzikadze on 3/31/17.
//

#include "config_based_getter.hpp"
#include "uniform_gene_chooser.hpp"
#include "custom_gene_chooser.hpp"


namespace ig_simulator {

AbstractVDJGeneChooserCPtr get_gene_chooser(const GeneChooserParams& config,
                                            const std::vector<germline_utils::CustomGeneDatabase>& db)
{
    if (config.method == GeneChooserMethod::Uniform)
        return AbstractVDJGeneChooserCPtr(new UniformVDJGeneChooser(db));
    // TODO add Custom
    // else if (config.method == GeneChooserMethod::Custom)
    //     return AbstractVDJGeneChooserCPtr(new CustomGeneChooser(db, config.custom_gene_chooser_params));
    VERIFY(false);
}

} // End namespace ig_simulator


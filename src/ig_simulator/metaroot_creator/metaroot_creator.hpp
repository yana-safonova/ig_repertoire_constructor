//
// Created by Andrew Bzikadze on 3/20/17.
//

#pragma once

#include "germline_utils/chain_type.hpp"
#include "gene_chooser/abstract_gene_chooser.hpp"
#include "nucleotides_remover/abstract_nucleotides_remover.hpp"
#include "p_nucleotides_creator/abstract_nucleotides_creator.hpp"
#include "n_nucleotides_inserter/abstract_n_nucleotides_inserter.hpp"

namespace ig_simulator {

class MetarootCreator {
private:
    const double prob_cleavage = 0.5;
    const bool is_vdj;
    AbstractVDJGeneChooserPtr gene_chooser_p;
    AbstractNucleotidesRemoverPtr nucl_remover_p;
    AbstractPNucleotidesCreatorPtr nucl_creator_p;
    AbstractNNucleotidesInserterPtr nucl_inserter_p;

public:
    MetarootCreator(AbstractVDJGeneChooserPtr gene_chooser_p_,
                    AbstractNucleotidesRemoverPtr nucl_remover_p_,
                    AbstractPNucleotidesCreatorPtr nucl_creator_p_,
                    AbstractNNucleotidesInserterPtr nucl_inserter_p_,
                    const germline_utils::ChainType& locus):
            is_vdj(locus.IsVDJ()),
            gene_chooser_p(std::move(gene_chooser_p_)),
            nucl_remover_p(std::move(nucl_remover_p_)),
            nucl_creator_p(std::move(nucl_creator_p_)),
            nucl_inserter_p(std::move(nucl_inserter_p_))
    { }

    void CreateRoot() const;
};

} // End namespace ig_simulator

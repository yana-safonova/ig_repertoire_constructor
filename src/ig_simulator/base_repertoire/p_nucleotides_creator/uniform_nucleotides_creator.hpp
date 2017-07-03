//
// Created by Andrew Bzikadze on 3/20/17.
//

#pragma once

#include "abstract_nucleotides_creator.hpp"
#include "ig_simulator_config.hpp"

namespace ig_simulator {

class UniformPNucleotidesCreator final : public AbstractPNucleotidesCreator {
private:
    size_t max_create_v_gene;
    size_t max_create_d_gene_left;
    size_t max_create_d_gene_right;
    size_t max_create_j_gene;

public:
    explicit UniformPNucleotidesCreator(
        const PNucleotidesCreatorParams::UniformCreatorParams config) :
            max_create_v_gene(config.max_create_v_gene),
            max_create_d_gene_left(config.max_create_d_gene_left),
            max_create_d_gene_right(config.max_create_d_gene_right),
            max_create_j_gene(config.max_create_j_gene)
    { }

    virtual size_t CreateInVGene()      const override;
    virtual size_t CreateInDGeneLeft()  const override;
    virtual size_t CreateInDGeneRight() const override;
    virtual size_t CreateInJGene()      const override;
};

} // End namespace ig_simulator

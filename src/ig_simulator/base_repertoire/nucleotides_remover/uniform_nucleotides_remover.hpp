//
// Created by Andrew Bzikadze on 3/17/17.
//

#pragma once

#include <random>
#include "abstract_nucleotides_remover.hpp"
#include "ig_simulator_config.hpp"

namespace ig_simulator {

class UniformNucleotidesRemover : public AbstractNucleotidesRemover {
private:
    size_t max_remove_v_gene;
    size_t max_remove_d_gene_left;
    size_t max_remove_d_gene_right;
    size_t max_remove_j_gene;

public:
    UniformNucleotidesRemover(
        const NucleotidesRemoverParams::UniformRemoverParams config) :
            AbstractNucleotidesRemover(),
            max_remove_v_gene(config.max_remove_v_gene),
            max_remove_d_gene_left(config.max_remove_d_gene_left),
            max_remove_d_gene_right(config.max_remove_d_gene_right),
            max_remove_j_gene(config.max_remove_j_gene)
    { }

    virtual size_t RemoveInVGene() const override;
    virtual size_t RemoveInDGeneLeft() const override;
    virtual size_t RemoveInDGeneRight() const override;
    virtual size_t RemoveInJGene() const override;

    virtual ~UniformNucleotidesRemover() { }
};

} // End namespace ig_simulator
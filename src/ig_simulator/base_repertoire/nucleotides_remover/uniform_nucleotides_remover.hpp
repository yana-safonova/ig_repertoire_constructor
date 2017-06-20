//
// Created by Andrew Bzikadze on 3/17/17.
//

#pragma once

#include <random>
#include "abstract_nucleotides_remover.hpp"
#include "ig_simulator_config.hpp"

namespace ig_simulator {

class UniformNucleotidesRemover final : public AbstractNucleotidesRemover {
private:
    const size_t max_remove_v_gene;
    const size_t max_remove_d_gene_left;
    const size_t max_remove_d_gene_right;
    const size_t max_remove_j_gene;

public:
    explicit UniformNucleotidesRemover(
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
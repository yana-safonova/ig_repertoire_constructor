//
// Created by Andrew Bzikadze on 3/17/17.
//

#pragma once

#include <random>
#include "abstract_nucleotides_remover.hpp"

namespace ig_simulator {

class UniformNucleotidesRemover : public AbstractNucleotidesRemover {
private:
    size_t max_remove_v_gene = 20;
    size_t max_remove_d_gene_left = 4;
    size_t max_remove_d_gene_right = 4;
    size_t max_remove_j_gene = 5;

public:
    UniformNucleotidesRemover() :
        AbstractNucleotidesRemover()
    { }

    virtual size_t RemoveInVGene() const override;
    virtual size_t RemoveInDGeneLeft() const override;
    virtual size_t RemoveInDGeneRight() const override;
    virtual size_t RemoveInJGene() const override;

    virtual ~UniformNucleotidesRemover() { }
};

} // End namespace ig_simulator
//
// Created by Andrew Bzikadze on 3/20/17.
//

#pragma once

#include "abstract_nucleotides_creator.hpp"

namespace ig_simulator {

class UniformPNucleotidesCreator : public AbstractPNucleotidesCreator {
private:
    size_t max_create_v_gene = 5;
    size_t max_create_d_gene_left = 3;
    size_t max_create_d_gene_right = 3;
    size_t max_create_j_gene = 3;

public:
    virtual size_t CreateInVGene()      const override;
    virtual size_t CreateInDGeneLeft()  const override;
    virtual size_t CreateInDGeneRight() const override;
    virtual size_t CreateInJGene()      const override;
};

} // End namespace ig_simulator

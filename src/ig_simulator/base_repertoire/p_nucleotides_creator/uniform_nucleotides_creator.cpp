//
// Created by Andrew Bzikadze on 3/20/17.
//

#include "uniform_nucleotides_creator.hpp"
#include "simulation_routines.hpp"

namespace ig_simulator {

size_t UniformPNucleotidesCreator::CreateInVGene()      const { return random_index(0, max_create_v_gene); }
size_t UniformPNucleotidesCreator::CreateInDGeneLeft()  const { return random_index(0, max_create_d_gene_left); }
size_t UniformPNucleotidesCreator::CreateInDGeneRight() const { return random_index(0, max_create_d_gene_right); }
size_t UniformPNucleotidesCreator::CreateInJGene()      const { return random_index(0, max_create_j_gene); }

} // End namespace ig_simulator

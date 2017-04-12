//
// Created by Andrew Bzikadze on 3/17/17.
//

#include "uniform_nucleotides_remover.hpp"
#include "simulation_routines.hpp"

namespace ig_simulator {

size_t UniformNucleotidesRemover::RemoveInVGene()      const { return random_index(0, max_remove_v_gene); }
size_t UniformNucleotidesRemover::RemoveInDGeneLeft()  const { return random_index(0, max_remove_d_gene_left); }
size_t UniformNucleotidesRemover::RemoveInDGeneRight() const { return random_index(0, max_remove_d_gene_right); }
size_t UniformNucleotidesRemover::RemoveInJGene()      const { return random_index(0, max_remove_j_gene); }

} // End namespace ig_simulator

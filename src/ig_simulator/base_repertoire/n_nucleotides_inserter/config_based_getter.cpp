//
// Created by Andrew Bzikadze on 3/31/17.
//

#include "config_based_getter.hpp"
#include "abstract_n_nucleotides_inserter.hpp"
#include "uniform_n_nucleotides_inserter.hpp"


namespace ig_simulator {

AbstractNNucleotidesInserterCPtr get_nucleotides_inserter(const NNucleotidesInserterParams & config)
{
    if (config.method == NNucleotidesInserterMethod::Uniform)
        return AbstractNNucleotidesInserterCPtr(new UniformNNucleotidesInserter(config.uniform_inserter_params));
    VERIFY(false);
}

} // End namespace ig_simulator


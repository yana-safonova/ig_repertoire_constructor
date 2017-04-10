//
// Created by Andrew Bzikadze on 3/31/17.
//

#include "config_based_getter.hpp"
#include "uniform_nucleotides_remover.hpp"


namespace ig_simulator {

AbstractNucleotidesRemoverCPtr get_nucleotides_remover(const NucleotidesRemoverParams & config)
{
    if (config.method == NucleotidesRemoverMethod::Uniform)
        return AbstractNucleotidesRemoverCPtr(new UniformNucleotidesRemover(config.uniform_remover_params));
    VERIFY(false);
}

} // End namespace ig_simulator


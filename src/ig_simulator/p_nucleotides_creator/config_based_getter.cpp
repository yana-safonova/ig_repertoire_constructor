//
// Created by Andrew Bzikadze on 3/31/17.
//

#include "config_based_getter.hpp"
#include "uniform_nucleotides_creator.hpp"


namespace ig_simulator {

AbstractPNucleotidesCreatorPtr get_nucleotides_creator(const PNucleotidesCreatorParams &config)
{
    if (config.method == PNucleotidesCreatorMethod::Uniform)
        return AbstractPNucleotidesCreatorPtr(new UniformPNucleotidesCreator(config.uniform_creator_params));
    VERIFY(false);
}

} // End namespace ig_simulator


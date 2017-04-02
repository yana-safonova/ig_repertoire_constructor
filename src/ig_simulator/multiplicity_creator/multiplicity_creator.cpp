//
// Created by Andrew Bzikadze on 3/27/17.
//

#include "multiplicity_creator.hpp"

namespace ig_simulator {

AbstractMultiplicityCreatorPtr get_multiplicity_creator(const MultiplicityCreatorParams &config) {
    if (config.method == MultiplicityCreatorMethod::Geometric) {
        return AbstractMultiplicityCreatorPtr(new GeometricMultiplicityCreator(config.geometric_params));
    }
    VERIFY(false);
}

} // End namespace ig_simulator

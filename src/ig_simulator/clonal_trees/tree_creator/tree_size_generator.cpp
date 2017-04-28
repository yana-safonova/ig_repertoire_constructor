//
// Created by Andrew Bzikadze on 4/11/17.
//

#include "tree_size_generator.hpp"

namespace ig_simulator {

AbstractTreeSizeGeneratorCPtr get_tree_size_generator(const TreeSizeGeneratorParams& config) {
    using Method = TreeSizeGeneratorParams::TreeSizeGeneratorMethod;
    if (config.method == Method::Geometric) {
        return AbstractTreeSizeGeneratorCPtr(new GeometricTreeSizeGenerator(config.geometric_params));
    }
    VERIFY(false);
}

} // End namespace ig_simulator

//
// Created by Andrew Bzikadze on 3/20/17.
//

#pragma once

#include <limits>
#include <random>
#include "random_generator.hpp"

namespace ig_simulator {

size_t random_index(size_t low = 0, size_t high = std::numeric_limits<size_t>::max());

template<typename FloatingPoint = double>
double uniform_floating_point(FloatingPoint low = 0., FloatingPoint high = 0.);

} // End namespace ig_simulator

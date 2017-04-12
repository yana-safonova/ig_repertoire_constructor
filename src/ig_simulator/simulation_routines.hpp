//
// Created by Andrew Bzikadze on 3/20/17.
//

#pragma once

#include <limits>
#include <random>
#include "random_generator.hpp"

namespace ig_simulator {

size_t random_index(size_t low = 0, size_t high = std::numeric_limits<size_t>::max());
double uniform_double(double low = 0., double high = 1.);

} // End namespace ig_simulator

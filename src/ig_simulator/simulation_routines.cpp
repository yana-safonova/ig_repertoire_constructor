//
// Created by Andrew Bzikadze on 3/20/17.
//

#include "simulation_routines.hpp"

namespace ig_simulator {

size_t random_index(size_t low, size_t high) {
    std::uniform_int_distribution<size_t> d(low, high);
    return d(MTSingleton::GetInstance());
}

template<typename FloatingPoint>
double uniform_floating_point(FloatingPoint low, FloatingPoint high) {
    static_assert(std::is_floating_point<FloatingPoint>::value, "Type has to be floating point");
    std::uniform_real_distribution<FloatingPoint> d(low, high);
    return d(MTSingleton::GetInstance());
}

} // End namespace ig_simulator

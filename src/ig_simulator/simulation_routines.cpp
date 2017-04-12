//
// Created by Andrew Bzikadze on 3/20/17.
//

#include "simulation_routines.hpp"

namespace ig_simulator {

size_t random_index(size_t low, size_t high) {
    std::uniform_int_distribution<size_t> d(low, high);
    return d(MTSingleton::GetInstance());
}

double uniform_double(double low, double high) {
    std::uniform_real_distribution<double> d(low, high);
    return d(MTSingleton::GetInstance());
}

} // End namespace ig_simulator

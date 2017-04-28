//
// Created by Andrew Bzikadze on 3/17/17.
//

#pragma once

#include <random>

namespace ig_simulator {

// This code is written after consulting with @eodus
template<class STLRandomGenerator, class Sseq=unsigned int>
class RandomGeneratorSingleton {
private:
    STLRandomGenerator generator_;

private:
    RandomGeneratorSingleton(Sseq seed=std::random_device()()) :
        generator_(seed)
    { }

public:
    static void SetSeed(Sseq seed = std::random_device()()) {
        RandomGeneratorSingleton::GetInstance().seed(seed);
    }

    static STLRandomGenerator& GetInstance() {
        static RandomGeneratorSingleton rg;
        return rg.generator_;
    }
};

using MTSingleton = RandomGeneratorSingleton<std::mt19937>;

} // End namespace ig_simulator

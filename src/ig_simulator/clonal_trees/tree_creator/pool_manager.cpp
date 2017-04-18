//
// Created by Andrew Bzikadze on 4/9/17.
//

#include "pool_manager.hpp"

namespace ig_simulator {

std::pair<size_t, bool> UniformPoolManager::GetIndex() {
    double raw_index = uniform_double(0, pool.Sum());
    size_t index, freq;
    std::tie(index, freq) = pool.LowerBound(raw_index);
    VERIFY(freq == 1);
    pool.Insert(max_index, 1);
    bool ret_to_pool = ret_to_pool_distr(MTSingleton::GetInstance());
    if (not ret_to_pool) {
        pool.Erase(index, freq);
    }
    max_index++;
    return { index, ret_to_pool };
}

std::pair<size_t, bool> WideTreePoolManager::GetIndex() {
    double raw_index = uniform_double(0, pool.Sum());
    size_t index, freq;
    std::tie(index, freq) = pool.LowerBound(raw_index);
    pool.Insert(max_index, 1);
    bool ret_to_pool = ret_to_pool_distr(MTSingleton::GetInstance());
    if (ret_to_pool) {
        pool.SetFreq(index, freq, freq + 1);
    } else {
        pool.Erase(index, freq);
    }
    max_index++;
    return { index, ret_to_pool };
}

std::pair<size_t, bool> DeepTreePoolManager::GetIndex() {
    double raw_index = uniform_double(0, pool.Sum());
    size_t index, freq;
    std::tie(index, freq) = pool.LowerBound(raw_index);
    pool.Insert(max_index, freq + 1);
    bool ret_to_pool = ret_to_pool_distr(MTSingleton::GetInstance());
    if (ret_to_pool) {
        pool.SetFreq(index, freq, freq + 1);
    } else {
        pool.Erase(index, freq);
    }
    max_index++;
    return { index, ret_to_pool };
}

} // End namespace ig_simulator

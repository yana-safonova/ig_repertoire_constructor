//
// Created by Andrew Bzikadze on 4/12/17.
//

#include "shm_creator.hpp"
#include "random_generator.hpp"

namespace ig_simulator {

Node::SHM_Vector PoissonShmCreator::GenerateSHM_Vector(size_t length) const {
    std::uniform_int_distribution<size_t> ind_distr(0, length);
    size_t mut_numb = distribution(MTSingleton::GetInstance()) + 1;
    std::vector<size_t> mut_inds;
    mut_inds.reserve(mut_numb);
    while(mut_inds.size() < mut_numb) {
        size_t ind = ind_distr(MTSingleton::GetInstance());
        if (std::find(mut_inds.begin(), mut_inds.end(), ind) == mut_inds.end()) {
            mut_inds.emplace_back(ind);
        }
    }

    std::uniform_int_distribution<size_t> mut_distr(0, 2);
    Node::SHM_Vector shm_vector;
    shm_vector.reserve(mut_numb);
    for(const auto& mut_ind : mut_inds) {
        shm_vector.emplace_back(mut_ind, mut_distr(MTSingleton::GetInstance()));
    }
    return shm_vector;
}

} // End namespace ig_simulator

//
// Created by Andrew Bzikadze on 4/12/17.
//

#include "shm_creator.hpp"
#include "random_generator.hpp"

namespace ig_simulator {

Node::SHM_Vector PoissonShmCreator::GenerateSHM_Vector(const std::string& seq) const {
    size_t length = seq.length();
    std::uniform_int_distribution<size_t> ind_distr(fix_left, length - 1 - fix_right);
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
        seqan::Dna5 old_nucl { seq[mut_ind] };
        size_t ind_nucl_old = old_nucl.value;
        size_t ind_nucl_new = mut_distr(MTSingleton::GetInstance());
        seqan::Dna new_nucl { ind_nucl_new < ind_nucl_old ? ind_nucl_new : ((ind_nucl_new + 1) & 3) };
        VERIFY_MSG(old_nucl != 'N' ? old_nucl != new_nucl : true,
                   std::string("Old nucl: ") << old_nucl
                       << ", New nucl: " << new_nucl
                       << " old nucl index: " << ind_nucl_old
                       << " new nucl index: " << ind_nucl_new
        );
        shm_vector.emplace_back(mut_ind, old_nucl, new_nucl);
    }
    return shm_vector;
}

AbstractShmCreatorCPtr get_shm_creator(const vj_finder::VJFinderConfig& vjf_config,
                                       const SHM_CreatorParams& config)
{
    using SHM_CreatorMethod = SHM_CreatorParams::SHM_CreatorMethod;
    if (config.method == SHM_CreatorMethod::Poisson) {
        return std::unique_ptr<AbstractShmCreator>(new PoissonShmCreator(vjf_config, config.poisson_params));
    }
    VERIFY(false);
}

} // End namespace ig_simulator

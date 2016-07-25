#include "recombination_estimator.hpp"
#include "../recombination/recombination.hpp"

class HcRecombinationEstimator : public RecombinationEstimator<HCRecombination> {
    std::vector<size_t> num_recombinations_;
    std::vector<size_t> min_number_shms_;
    std::vector<int> v_event_lens_;
    std::vector<int> d_levent_lens_;
    std::vector<int> d_revent_lens_;
    std::vector<int> j_event_lens_;

public:
    void Update(HcRecombinationStoragePtr recombination_storage);

    void OutputRecombinationNumber();

    void OutputSHMsDistribution();

    void OutputRecombinationEvents();
};
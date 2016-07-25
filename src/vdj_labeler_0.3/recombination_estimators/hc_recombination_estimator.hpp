#include "recombination_estimator.hpp"
#include "recombination_utils/recombination.hpp"

namespace vdj_labeler {

class HcRecombinationEstimator: public RecombinationEstimator<recombination_utils::HCRecombination> {
    std::vector<size_t> num_recombinations_;
    std::vector<size_t> min_number_shms_;
    std::vector<int> v_event_lens_;
    std::vector<int> d_levent_lens_;
    std::vector<int> d_revent_lens_;
    std::vector<int> j_event_lens_;

public:
    void Update(const recombination_utils::HcRecombinationStorage &recombination_storage) override;

    void OutputRecombinationNumber() const override;

    void OutputSHMsDistribution() const override;

    void OutputRecombinationEvents() const;
};

} // End namespace vdj_labeler
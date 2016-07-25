#include "recombination_utils/recombination_storage.hpp"

namespace vdj_labeler {

template<class Recombination>
class RecombinationEstimator {
public:
    RecombinationEstimator() = default;

    RecombinationEstimator(const RecombinationEstimator &) = delete;
    RecombinationEstimator &operator=(const RecombinationEstimator &) = delete;
    RecombinationEstimator(RecombinationEstimator &&) = delete;
    RecombinationEstimator &operator=(RecombinationEstimator &&)      = delete;

    virtual void Update(const recombination_utils::HcRecombinationStorage &recombination_storage) = 0;

    virtual void OutputRecombinationNumber() const = 0;

    virtual void OutputSHMsDistribution() const = 0;
};

} // End namespace vdj_labeler
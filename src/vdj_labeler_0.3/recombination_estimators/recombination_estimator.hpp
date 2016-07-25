#include "../recombination_generation/recombination_storage.hpp"

template<class Recombination>
class RecombinationEstimator {
public:
    virtual void Update(std::shared_ptr<RecombinationStorage<Recombination> > recombination_storage) = 0;

    virtual void OutputRecombinationNumber() = 0;

    virtual void OutputSHMsDistribution() = 0;
};
#include "clusterer.hpp"

namespace clusterer {

    const ClusteringMode ClusteringMode::hamming = ClusteringMode(
            [](const seqan::Dna5String &first, const seqan::Dna5String &second) {
                return static_cast<size_t>(-half_sw_banded(first, second, 0, -1, -1, [](int) -> int { return 0; }, 0));
            }, 10);

    const ClusteringMode ClusteringMode::edit = ClusteringMode(
            [](const seqan::Dna5String &first, const seqan::Dna5String &second) {
                return static_cast<unsigned long>(get_sw_dist(first, second));
            }, 10);


    ReflexiveUmiPairsIterator ReflexiveUmiPairsIterator::operator++() {
        current_ ++;
        return *this;
    }

    ReflexiveUmiPairsIterator ReflexiveUmiPairsIterator::operator++(int) {
        const ReflexiveUmiPairsIterator itr(current_, last_);
        current_ ++;
        return itr;
    }

    std::pair<size_t, size_t> ReflexiveUmiPairsIterator::operator*() const {
        return std::make_pair(current_, current_);
    }

    bool ReflexiveUmiPairsIterator::operator==(ReflexiveUmiPairsIterator other) const {
        VERIFY_MSG(other.last_ == last_, "Comparing different iterators.");
        return current_ == other.current_;
    }

    bool ReflexiveUmiPairsIterator::operator!=(ReflexiveUmiPairsIterator other) const {
        return !(*this == other);
    }

    ReflexiveUmiPairsIterator ReflexiveUmiPairsIterable::begin() const {
        return ReflexiveUmiPairsIterator(0, count_);
    }

    ReflexiveUmiPairsIterator ReflexiveUmiPairsIterable::end() const {
        return ReflexiveUmiPairsIterator(count_, count_);
    }
}

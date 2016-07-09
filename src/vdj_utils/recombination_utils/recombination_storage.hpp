#pragma once

#include "recombination.hpp"
#include "read_archive.hpp"

#include <vector>

namespace recombination_utils {

template<class Recombination>
class RecombinationStorage {
    core::ReadPtr read_ptr_;
    std::vector <Recombination> recombinations_;

    bool CheckConsistency(Recombination recombination) {
        return recombination.ReadId() == read_ptr_->id;
    }

public:
    RecombinationStorage() = delete;

    explicit RecombinationStorage(const core::ReadPtr read_ptr) :
        read_ptr_(read_ptr) { }

    explicit RecombinationStorage(const RecombinationStorage&) = default;

    void AddRecombination(const Recombination &recombination) {
        if (CheckConsistency(recombination))
            recombinations_.push_back(recombination);
    }

    typedef typename std::vector<Recombination>::const_iterator recombination_iterator;

    recombination_iterator cbegin() const { return recombinations_.cbegin(); }

    recombination_iterator cend() const { return recombinations_.cend(); }

    size_t size() const { return recombinations_.size(); }

    core::ReadPtr Read() const { return read_ptr_; }

    const Recombination& operator[](size_t index) const {
        assert(index < size());
        return recombinations_[index];
    }

    Recombination& operator[](size_t index) {
        assert(index < size());
        return recombinations_[index];
    }

    size_t ReadId() const {
        assert(read_ptr != nullptr);
        return read_ptr_-> id;
    }
};

typedef RecombinationStorage<HCRecombination> HcRecombinationStorage;

typedef std::shared_ptr <HcRecombinationStorage> HcRecombinationStoragePtr;

} // End namespace recombination_utils
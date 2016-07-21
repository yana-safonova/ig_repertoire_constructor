#pragma once

#include "recombination.hpp"
#include "read_archive.hpp"

#include <vector>

namespace recombination_utils {

template<class Recombination>
class RecombinationStorage {
    const core::Read* read_ptr_;
    std::vector <Recombination> recombinations_;

    bool CheckConsistency(const Recombination &recombination) {
        return recombination.ReadId() == read_ptr_->id;
    }

public:
    RecombinationStorage() :
        read_ptr_(nullptr),
        recombinations_()
    { }

    explicit RecombinationStorage(const core::Read* read_ptr) :
        read_ptr_(read_ptr)
    { }

    RecombinationStorage(const RecombinationStorage&) = default;
    RecombinationStorage(RecombinationStorage&&) = default;
    RecombinationStorage& operator=(const RecombinationStorage&) = default;
    RecombinationStorage& operator=(RecombinationStorage&&) = default;

    void AddRecombination(Recombination recombination) {
        if (CheckConsistency(recombination))
            recombinations_.emplace_back(std::move(recombination));
    }

    typedef typename std::vector<Recombination>::iterator recombination_iterator;
    typedef typename std::vector<Recombination>::const_iterator recombination_const_iterator;

    recombination_iterator       begin ()       { return recombinations_.begin (); }
    recombination_const_iterator begin () const { return recombinations_.begin (); }
    recombination_const_iterator cbegin() const { return recombinations_.cbegin(); }
    recombination_iterator       end   ()       { return recombinations_.end   (); }
    recombination_const_iterator end   () const { return recombinations_.end   (); }
    recombination_const_iterator cend  () const { return recombinations_.cend  (); }

    size_t size() const { return recombinations_.size(); }

    const core::Read* ReadPtr() const { return read_ptr_; }
    const core::Read& Read() const {
        assert(read_ptr_ != nullptr);
        return *read_ptr_;
    }

    // Andrey: Users should not be allowed to change the elements. Thus, the other version is ommited:
    // Recombination& operator[](size_t index) {
    //     assert(index < size());
    //     return recombinations_[index];
    // }
    const Recombination& operator[](size_t index) const {
        assert(index < size());
        return recombinations_[index];
    }

    size_t ReadId() const {
        assert(read_ptr_ != nullptr);
        return read_ptr_-> id;
    }
};

typedef RecombinationStorage<HCRecombination> HcRecombinationStorage;

typedef std::shared_ptr <HcRecombinationStorage> HcRecombinationStoragePtr;

} // End namespace recombination_utils
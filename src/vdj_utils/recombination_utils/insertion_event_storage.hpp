#pragma once

#include "nongenomic_insertion.hpp"

#include <vector>

namespace recombination_utils {

class InsertionEventStorage {
    std::vector<NongenomicInsertion> insertions_;

public:
    InsertionEventStorage() = default;
    InsertionEventStorage(const InsertionEventStorage&) = default;
    InsertionEventStorage(InsertionEventStorage&&) = default;
    InsertionEventStorage& operator=(const InsertionEventStorage&) = default;
    InsertionEventStorage& operator=(InsertionEventStorage&&) = default;

    explicit InsertionEventStorage(std::vector<NongenomicInsertion> insertions) :
        insertions_(std::move(insertions))
    { }

    void AddInsertion(NongenomicInsertion insertion) {
        insertions_.emplace_back(std::move(insertion));
    }

    typedef std::vector<NongenomicInsertion>::iterator insertion_event_iterator;
    typedef std::vector<NongenomicInsertion>::const_iterator insertion_event_const_iterator;

    insertion_event_iterator       begin ()       { return insertions_.begin();  }
    insertion_event_const_iterator begin () const { return insertions_.begin();  }
    insertion_event_const_iterator cbegin() const { return insertions_.cbegin(); }
    insertion_event_iterator       end   ()       { return insertions_.end  ();  }
    insertion_event_const_iterator end   () const { return insertions_.end  ();  }
    insertion_event_const_iterator cend  () const { return insertions_.cend  (); }

    size_t size() const { return insertions_.size(); }
};

typedef std::shared_ptr<InsertionEventStorage> InsertionEventStoragePtr;

} // End namespace recombination_utils

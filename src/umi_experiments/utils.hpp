#pragma once

#include <unordered_set>
#include <logger/log_writers.hpp>
#include <seqan/seq_io.h>

void create_console_logger();

size_t get_sw_dist(const seqan::Dna5String& first, const seqan::Dna5String& second);

template <typename From, typename To>
class ManyToManyCorrespondence {
public:
    ManyToManyCorrespondence() = default;
    ManyToManyCorrespondence(const ManyToManyCorrespondence& other);

    const std::unordered_set<To>& forth(const From& from) const { return forth_.find(from)->second; }
    const std::unordered_set<From>& back(const To& to) const { return back_.find(to)->second; }

    // returns an element from the 'to' set which is equal to 'to' parameter
    const To getTo(const To& to) const { return back_.find(to)->first; }
    size_t toSize() const { return back_.size(); }
    const std::unordered_set<To> toSet() const;

    // returns true if the 'to' parameter is presented by some links
    bool removeTo(const To &to);
    void add(const From& from, const To& to);
    void add(const std::unordered_set<From>& from_set, const To& to);

private:
    std::unordered_map<From, std::unordered_set<To>> forth_;
    std::unordered_map<To, std::unordered_set<From>> back_;
};

template <typename From, typename To>
ManyToManyCorrespondence<From, To>::ManyToManyCorrespondence(const ManyToManyCorrespondence& other) : forth_(other.forth_), back_(other.back_) {
    VERIFY_MSG(forth_.size() == other.forth_.size() && back_.size() == other.back_.size(), "!!!");
}

template <typename From, typename To>
const std::unordered_set<To> ManyToManyCorrespondence<From, To>::toSet() const {
    std::unordered_set<To> result;
    for (const auto& entry : back_) {
        result.insert(entry.first);
    }
    return result;
};

template <typename From, typename To>
bool ManyToManyCorrespondence<From, To>::removeTo(const To &to) {
    size_t contains = back_.count(to);
    if (contains == 0) return false;
    for (const auto& from : back_[to]) {
        forth_[from].erase(to);
    }
    back_.erase(to);
    return true;
}

template <typename From, typename To>
void ManyToManyCorrespondence<From, To>::add(const From& from, const To& to) {
    VERIFY_MSG(forth_[from].insert(to).second, "Adding already existing mapping.");
    VERIFY_MSG(back_[to].insert(from).second, "Adding already existing mapping.");
}

template <typename From, typename To>
void ManyToManyCorrespondence<From, To>::add(const std::unordered_set<From>& from_set, const To& to) {
    for (const auto& from : from_set) {
        add(from, to);
    }
}

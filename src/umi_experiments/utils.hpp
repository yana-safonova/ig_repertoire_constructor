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

    size_t rightSize() const { return back_.size(); }

    // returns true if the to parameter is presented by some links
    bool removeSecond(const To& to);
    void add(const From& from, const To& to);
    void add(const std::unordered_set<From>& from_set, const To& to);

private:
    std::unordered_map<From, std::unordered_set<To>> forth_;
    std::unordered_map<To, std::unordered_set<From>> back_;
};

template <typename From, typename To>
ManyToManyCorrespondence<From, To>::ManyToManyCorrespondence(const ManyToManyCorrespondence& other) : forth_(other.forth_), back_(other.back_) {}

template <typename From, typename To>
bool ManyToManyCorrespondence<From, To>::removeSecond(const To& to) {
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
    forth_[from].insert(to);
    back_[to].insert(from);
}

template <typename From, typename To>
void ManyToManyCorrespondence<From, To>::add(const std::unordered_set<From>& from_set, const To& to) {
    for (const auto& from : from_set) {
        add(from, to);
    }
}

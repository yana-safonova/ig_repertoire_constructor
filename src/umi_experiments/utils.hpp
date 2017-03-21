#pragma once

#include <unordered_set>
#include <logger/log_writers.hpp>
#include <seqan/seq_io.h>
#include "../ig_tools/utils/string_tools.hpp"

void create_console_logger();

size_t get_sw_dist(const seqan::Dna5String& first, const seqan::Dna5String& second);

template <typename From, typename To, typename FromHash, typename FromEquals, typename ToHash, typename ToEquals>
class ManyToManyCorrespondence {
public:
    ManyToManyCorrespondence() = default;
    ManyToManyCorrespondence(const FromHash& from_hash, const FromEquals& from_equals, const ToHash& to_hash, const ToEquals& to_equals);
    ManyToManyCorrespondence(const ManyToManyCorrespondence& other);

    const std::unordered_set<To, ToHash, ToEquals>& forth(const From& from) const {
        VERIFY(forth_.find(from) != forth_.end());
        return forth_.find(from)->second;
    }
    const std::unordered_set<From, FromHash, FromEquals>& back(const To& to) const {
        VERIFY(back_.find(to) != back_.end());
        return back_.find(to)->second;
    }

    // returns an element from the 'to' set which is equal to 'to' parameter
    const To getTo(const To& to) const {
        const auto& itr = back_.find(to);
        VERIFY_MSG(itr != back_.end(), "Not found.");
        return itr->first;
    }
    size_t toSize() const { return back_.size(); }
    const std::unordered_set<To, ToHash, ToEquals> toSet() const;

    const std::unordered_set<From, FromHash, FromEquals> fromSet() const;

    // returns true if the 'to' parameter is presented by some links
    bool removeTo(const To &to);
    void add(const From& from, const To& to);
    void add(const std::unordered_set<From, FromHash, FromEquals>& from_set, const To& to);

private:
    std::unordered_map<From, std::unordered_set<To, ToHash, ToEquals>, FromHash, FromEquals> forth_;
    std::unordered_map<To, std::unordered_set<From, FromHash, FromEquals>, ToHash, ToEquals> back_;

    void check_state();
    DECL_LOGGER("ManyToManyCorrespondence");
};

template <typename From, typename To, typename FromHash, typename FromEquals, typename ToHash, typename ToEquals>
ManyToManyCorrespondence<From, To, FromHash, FromEquals, ToHash, ToEquals>::ManyToManyCorrespondence
        (const FromHash& from_hash, const FromEquals& from_equals, const ToHash& to_hash, const ToEquals& to_equals)
        : forth_(std::unordered_map<From, std::unordered_set<To, ToHash, ToEquals>, FromHash, FromEquals>(1, from_hash, from_equals)),
          back_(std::unordered_map<To, std::unordered_set<From, FromHash, FromEquals>, ToHash, ToEquals>(1, to_hash, to_equals)) {
    check_state();
}

template <typename From, typename To, typename FromHash, typename FromEquals, typename ToHash, typename ToEquals>
ManyToManyCorrespondence<From, To, FromHash, FromEquals, ToHash, ToEquals>::ManyToManyCorrespondence
        (const ManyToManyCorrespondence& other) : forth_(other.forth_), back_(other.back_) {
    check_state();
}

template <typename From, typename To, typename FromHash, typename FromEquals, typename ToHash, typename ToEquals>
const std::unordered_set<To, ToHash, ToEquals> ManyToManyCorrespondence<From, To, FromHash, FromEquals, ToHash, ToEquals>::toSet() const {
    std::unordered_set<To, ToHash, ToEquals> result;
    for (const auto& entry : back_) {
        result.insert(entry.first);
    }
    return result;
};

template <typename From, typename To, typename FromHash, typename FromEquals, typename ToHash, typename ToEquals>
const std::unordered_set<From, FromHash, FromEquals> ManyToManyCorrespondence<From, To, FromHash, FromEquals, ToHash, ToEquals>::fromSet() const {
    std::unordered_set<From, FromHash, FromEquals> result;
    for (const auto& entry : forth_) {
        result.insert(entry.first);
    }
    return result;
};

template <typename From, typename To, typename FromHash, typename FromEquals, typename ToHash, typename ToEquals>
bool ManyToManyCorrespondence<From, To, FromHash, FromEquals, ToHash, ToEquals>::removeTo(const To &to) {
    size_t contains = back_.count(to);
    if (contains == 0) {
        return false;
    }
    for (const auto& from : back_[to]) {
        forth_[from].erase(to);
        // should not be needed for clusterer, but could be expected in general
        if (forth_[from].empty()) {
            forth_.erase(from);
        }
    }
    back_.erase(to);
    VERIFY(back_.count(to) == 0);

//    check_state();
    return true;
}

template <typename From, typename To, typename FromHash, typename FromEquals, typename ToHash, typename ToEquals>
void ManyToManyCorrespondence<From, To, FromHash, FromEquals, ToHash, ToEquals>::add(const From& from, const To& to) {
    VERIFY_MSG(forth_[from].insert(to).second, "Adding already existing mapping.");
    VERIFY_MSG(back_[to].insert(from).second, "Adding already existing mapping.");
//    check_state();
}

template <typename From, typename To, typename FromHash, typename FromEquals, typename ToHash, typename ToEquals>
void ManyToManyCorrespondence<From, To, FromHash, FromEquals, ToHash, ToEquals>::add(const std::unordered_set<From, FromHash, FromEquals>& from_set, const To& to) {
    for (const auto& from : from_set) {
        add(from, to);
    }
    VERIFY_MSG(back_.find(to) != back_.end(), "Not added actually");
}

template <typename From, typename To, typename FromHash, typename FromEquals, typename ToHash, typename ToEquals>
void ManyToManyCorrespondence<From, To, FromHash, FromEquals, ToHash, ToEquals>::check_state() {
    {
        std::unordered_set<To, ToHash, ToEquals> forth_targets;
        for (const auto& entry : forth_) {
            forth_targets.insert(entry.second.begin(), entry.second.end());
        }
        VERIFY_MSG(forth_targets.size() == back_.size(), "To set of size " << back_.size() << ", but forth map points to " << forth_targets.size() << " targets");
    }
    {
        std::unordered_set<From, FromHash, FromEquals> back_targets;
        for (const auto& entry : back_) {
            back_targets.insert(entry.second.begin(), entry.second.end());
        }
        VERIFY_MSG(back_targets.size() == forth_.size(), "From set of size " << forth_.size() << ", but back map points to " << back_targets.size() << " targets");
    }

    for (const auto& entry : forth_) {
        for (const auto& to : entry.second) {
            VERIFY_MSG(back_.count(to), seqan_string_to_string(entry.first->GetString()) << " umi maps to cluster with center "
                    << seqan_string_to_string(to->GetSequence()) << ", but there's no such cluster in back map");
        }
    }
    for (const auto& entry : back_) {
        for (const auto& from : entry.second) {
            VERIFY_MSG(forth_.count(from), seqan_string_to_string(entry.first->GetSequence()) << "-centered cluster maps to umi "
                    << seqan_string_to_string(from->GetString()) << ", but there's no such umi in forth map");
        }
    }
}

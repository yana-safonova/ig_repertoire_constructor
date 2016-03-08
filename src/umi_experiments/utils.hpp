#pragma once

#include <unordered_set>
#include <logger/log_writers.hpp>
#include <seqan/seq_io.h>

void create_console_logger();

size_t get_sw_dist(const seqan::Dna5String& first, const seqan::Dna5String& second);

template <typename From, typename To>
class ManyToManyCorrespondence {
private:
    typedef std::shared_ptr<From> FromPtr;
    typedef std::shared_ptr<To> ToPtr;
public:
    ManyToManyCorrespondence() = default;
    ManyToManyCorrespondence(const ManyToManyCorrespondence& other);

    std::unordered_set<ToPtr> forth(FromPtr from) const { return forth_[from]; }
    std::unordered_set<FromPtr> back(ToPtr to) const { return back_[to]; }

    // returns true if the to parameter is presented by some links
    bool removeSecond(ToPtr to);
    void add(const FromPtr& from, const ToPtr& to);
    void add(const std::unordered_set<FromPtr>& from_set, const ToPtr& to);

private:
    std::unordered_map<FromPtr, std::unordered_set<ToPtr>> forth_;
    std::unordered_map<ToPtr, std::unordered_set<FromPtr>> back_;
};
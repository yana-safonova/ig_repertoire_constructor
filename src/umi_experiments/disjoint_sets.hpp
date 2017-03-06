#pragma once

#include <unordered_map>
#include <verify.hpp>

template <typename T, typename Hash = std::hash<T>>
class DisjointSets {
public:

    void addNewSet(const T element) {
        VERIFY_MSG(!parent_.count(element), "There's already such element somewhere in the sets.");
        parent_[element] = element;
        depth_[element] = 0;
    }

    T findRoot(const T& element) {
        if (parent_[element] == element) {
            return element;
        }
        return parent_[element] = findRoot(parent_[element]);
    }

    bool unite(const T& first, const T& second) {
        T first_parent = findRoot(first);
        T second_parent = findRoot(second);
        if (first_parent == second_parent) {
            return false;
        }
        if (depth_[first_parent] > depth_[second_parent]) {
            parent_[second_parent] = first_parent;
        } else {
            parent_[first_parent] = second_parent;
            if (depth_[first_parent] == depth_[second_parent]) {
                depth_[second_parent] ++;
            }
        }
        return true;
    }

private:
    std::unordered_map<T, T, Hash> parent_;
    std::unordered_map<T, size_t, Hash> depth_;
};
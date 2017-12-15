#pragma once

#include <logger/logger.hpp>

class Permutation {
    size_t num_vertices_;
    std::vector<size_t> direct_;
    std::vector<size_t> reverse_;

public:
    Permutation(size_t num_vertices) :
        num_vertices_(num_vertices) {
        for(size_t i = 0; i < num_vertices_; i++) {
            direct_.push_back(i);
            reverse_.push_back(i);
        }
    }

    void ReadFromFile(std::string filename);

    size_t Size() const { return num_vertices_; }

    const std::vector<size_t>& Direct() const { return direct_; }

    const std::vector<size_t>& Reverse() const { return reverse_; }

private:
    DECL_LOGGER("Permutation");
};

std::ostream& operator<<(std::ostream &out, const Permutation &permutation);

typedef std::shared_ptr<Permutation> PermutationPtr;

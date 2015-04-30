#pragma once

#include "../../include_me.hpp"

class Permutation {
    size_t num_vertices_;
    vector<size_t> direct_;
    vector<size_t> reverse_;

public:
    Permutation(size_t num_vertices) :
        num_vertices_(num_vertices) {
        for(size_t i = 0; i < num_vertices_; i++) {
            direct_.push_back(i);
            reverse_.push_back(i);
        }
    }

    void ReadFromFile(string filename) {
        ifstream perm_fhandler(filename.c_str());
        size_t index1 = 0;
        size_t index2 = 0;
        while(!perm_fhandler.eof()) {
            string line;
            getline(perm_fhandler, line);
            if(line == "")
                break;

            stringstream ss(line);
            ss >> index2;

            direct_[index1] = index2;
            reverse_[index2] = index1;
            index1++;
        }
        perm_fhandler.close();
    }

    size_t Size() const { return num_vertices_; }

    const vector<size_t>& Direct() const { return direct_; }

    const vector<size_t>& Reverse() const { return reverse_; }

private:
    DECL_LOGGER("Permutation");
};

ostream& operator<<(ostream &out, const Permutation &permutation) {
    out << "Direct: ";
    for(size_t i = 0; i < permutation.Size(); i++)
        out << i << ": " << permutation.Direct()[i] << ", ";
    out << endl;
    out << "Reverse: ";
    for(size_t i = 0; i < permutation.Size(); i++)
        out << i << ": " << permutation.Reverse()[i] << ", ";
    return out;
}

typedef shared_ptr<Permutation> PermutationPtr;
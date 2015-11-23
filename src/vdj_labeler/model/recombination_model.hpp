#pragma once

#include "standard.hpp"

#include <map>
#include <string>
#include <fstream>


using std::string;
using std::map;
using std::vector;

class IgGeneProbabilityModel {
    map<string, double> ig_gene_probabilities_;

public:
    IgGeneProbabilityModel() = delete;

    IgGeneProbabilityModel(map<string, double>);

    IgGeneProbabilityModel(const IgGeneProbabilityModel&) = default;

    IgGeneProbabilityModel(IgGeneProbabilityModel&&) = default;

    IgGeneProbabilityModel& operator=(const IgGeneProbabilityModel&) = default; 

    IgGeneProbabilityModel& operator=(IgGeneProbabilityModel&&) = default;

    virtual ~IgGeneProbabilityModel() = default;

    map<string, double>& GetGeneProbabilities() { return ig_gene_probabilities_; }
    const map<string, double>& GetGeneProbabilities() const {
        return ig_gene_probabilities_; 
    }

    void SetGeneProbabilities(map<string, double>& ig_gene_probabilities) {
        ig_gene_probabilities_ = ig_gene_probabilities;
    }

    // double operator[](const string& index_string) { return ig_gene_probabilities_[index_string]; }
    double GetProbabilityByGeneName(const string&) const;

    size_t size() const;

    IgGeneProbabilityModel(std::ifstream&);

    void print(std::ostream& out);
};



class NongenomicInsertionModel {
    static const size_t alphabet_size_;
    vector<double> insertion_probabilities_;
    vector<vector<double>> transition_matrix_; // (4, vector<double>(alphabet_size));

public:
    NongenomicInsertionModel() = delete;

    NongenomicInsertionModel(vector<double>, vector<vector<double>>);

    NongenomicInsertionModel(const NongenomicInsertionModel&) = default;

    NongenomicInsertionModel(NongenomicInsertionModel&&) = default;

    NongenomicInsertionModel& operator=(const NongenomicInsertionModel&) = default;

    NongenomicInsertionModel& operator=(NongenomicInsertionModel&&) = default;

    virtual ~NongenomicInsertionModel() = default;

    vector<double> GetInsertionProbabilities() { return insertion_probabilities_; }
    const vector<double>& GetInsertionProbabilities() const { return insertion_probabilities_; }

    vector<vector<double>> GetTransitionMatrix() { return transition_matrix_; }
    const vector<vector<double>>& GetTransitionMatrix() const { return transition_matrix_; }

    void SetInsertionProbabilities(const vector<double>& insertion_probabilities) {
        insertion_probabilities_ = insertion_probabilities;
    }

    void SetTransitionMatrix(const vector<vector<double>>& transition_matrix) {
        transition_matrix_ = transition_matrix;
    }

    NongenomicInsertionModel(std::ifstream&);

    void print(std::ostream& out);
};

class PalindromeDeletionModel {
    using DeletionTableMap = map<string, map<int, double>>;
    DeletionTableMap deletion_table_;

public:
    PalindromeDeletionModel() = delete;

    PalindromeDeletionModel(DeletionTableMap deletion_table);

    PalindromeDeletionModel(const PalindromeDeletionModel&) = default;

    PalindromeDeletionModel(PalindromeDeletionModel&&) = default;

    PalindromeDeletionModel& operator=(const PalindromeDeletionModel&) = default;

    PalindromeDeletionModel& operator=(PalindromeDeletionModel&&) = default;

    virtual ~PalindromeDeletionModel() = default;

    DeletionTableMap GetDeletionTable() { return deletion_table_; }
    const DeletionTableMap& GetDeletionTable() const { return deletion_table_; }

    void SetDeletionTable(const DeletionTableMap& deletion_table) {
        deletion_table_ = deletion_table;
    }

    PalindromeDeletionModel(std::ifstream&);

    void print(std::ostream& out);

};

class HCProbabilityRecombinationModel {
    IgGeneProbabilityModel V_gene_probability_model_;
    IgGeneProbabilityModel D_gene_probability_model_;
    IgGeneProbabilityModel J_gene_probability_model_;

    NongenomicInsertionModel VD_nongenomic_insertion_model_;
    NongenomicInsertionModel DJ_nongenomic_insertion_model_;

    PalindromeDeletionModel V_palindrome_deletion_model_;
    PalindromeDeletionModel J_palindrome_deletion_model_;
    PalindromeDeletionModel DLeft_palindrome_deletion_model_;
    PalindromeDeletionModel DRight_palindrome_deletion_model_;

public:
    HCProbabilityRecombinationModel() = delete;

    HCProbabilityRecombinationModel(std::ifstream&);

    void print(std::ostream& out);
};

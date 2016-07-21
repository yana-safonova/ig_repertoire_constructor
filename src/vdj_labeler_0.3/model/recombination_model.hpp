#pragma once

#include "germline_utils/germline_databases/immune_gene_database.hpp"
#include "germline_utils/germline_databases/chain_database.hpp"
#include "germline_utils/germline_gene_type.hpp"

#include <fstream>
#include <vector>
#include <utility>
#include <memory>

namespace vdj_labeler {

// Here we use only Heavy Chain (HC).
class IgGeneProbabilityModel {
    using IgGeneProbabilityVector = std::vector<double>;
    IgGeneProbabilityVector ig_gene_probabilities_;
    const germline_utils::ImmuneGeneDatabase* ig_gene_database_;

public:
    IgGeneProbabilityModel() = delete;

    IgGeneProbabilityModel(const IgGeneProbabilityVector&, const germline_utils::ImmuneGeneDatabase*);

    IgGeneProbabilityModel(const IgGeneProbabilityModel &) = default;
    IgGeneProbabilityModel(IgGeneProbabilityModel &&) = default;
    IgGeneProbabilityModel &operator=(const IgGeneProbabilityModel &) = default;
    IgGeneProbabilityModel &operator=(IgGeneProbabilityModel &&) = default;

    virtual ~IgGeneProbabilityModel() = default;

    const IgGeneProbabilityVector &GetIgGeneProbabilities() const { return ig_gene_probabilities_; }

    const germline_utils::ImmuneGeneDatabase* GetIgGeneDatabasePtr() const { return ig_gene_database_; }

    const germline_utils::ImmuneGeneDatabase &GetIgGeneDatabase() const {
        assert(ig_gene_database_ != nullptr);
        return *ig_gene_database_;
    }

    void SetGeneProbabilities(const IgGeneProbabilityVector &ig_gene_probabilities) {
        ig_gene_probabilities_ = ig_gene_probabilities;
    }

    using iterator = IgGeneProbabilityVector::iterator;
    using citerator = IgGeneProbabilityVector::const_iterator;

    iterator  begin ()       { return ig_gene_probabilities_.begin(); }
    citerator begin () const { return ig_gene_probabilities_.begin(); }
    citerator cbegin() const { return ig_gene_probabilities_.cbegin(); }
    iterator  end   ()       { return ig_gene_probabilities_.end(); }
    citerator end   () const { return ig_gene_probabilities_.end(); }
    citerator cend  () const { return ig_gene_probabilities_.cend(); }

    size_t size() const;

    IgGeneProbabilityModel(std::ifstream &, const germline_utils::ImmuneGeneDatabase*);

    IgGeneProbabilityModel(std::ifstream &, const germline_utils::ImmuneGeneDatabase &);

    //double GetProbabilityByGenId(const IgGene&) const;
    double GetProbabilityByGenId(const size_t &) const;
};

std::ostream &operator<<(std::ostream &, const IgGeneProbabilityModel &);

/***************************************************************************************************/

class NongenomicInsertionModel {
    std::vector<double> insertion_probabilities_;
    using NongenomicInsertionMatrix = std::vector<std::vector<double>>;
    NongenomicInsertionMatrix transition_matrix_;

public:
    NongenomicInsertionModel() = delete;

    NongenomicInsertionModel(std::vector<double>, NongenomicInsertionMatrix);

    NongenomicInsertionModel(const NongenomicInsertionModel &) = default;

    NongenomicInsertionModel(NongenomicInsertionModel &&) = default;

    NongenomicInsertionModel &operator=(const NongenomicInsertionModel &) = default;

    NongenomicInsertionModel &operator=(NongenomicInsertionModel &&) = default;

    virtual ~NongenomicInsertionModel() = default;

    const std::vector<double> &GetInsertionProbabilities() const { return insertion_probabilities_; }

    double GetInsertionProbabilityByLength(const long unsigned int) const;

    const NongenomicInsertionMatrix &GetTransitionMatrix() const { return transition_matrix_; }

    void SetInsertionProbabilities(const std::vector<double> &insertion_probabilities) {
        insertion_probabilities_ = insertion_probabilities;
    }

    void SetTransitionMatrix(const NongenomicInsertionMatrix &transition_matrix) {
        transition_matrix_ = transition_matrix;
    }

    explicit NongenomicInsertionModel(std::ifstream &);

    double GetTransitionProbability(char, char) const;

    double GetTransitionProbability(const std::pair<char, char> &) const;
};

std::ostream &operator<<(std::ostream &, const NongenomicInsertionModel &);

/***************************************************************************************************/

class PalindromeDeletionModel {
    using DeletionTableVector = std::vector<std::vector<double>>;
    DeletionTableVector deletion_table_;
    std::vector<int> deletion_length_;
    const germline_utils::ImmuneGeneDatabase* ig_gene_database_;

public:
    PalindromeDeletionModel() = delete;

    PalindromeDeletionModel(const DeletionTableVector &,
                            const std::vector<int> &,
                            const germline_utils::ImmuneGeneDatabase*);

    PalindromeDeletionModel(const PalindromeDeletionModel &) = default;

    PalindromeDeletionModel(PalindromeDeletionModel &&) = default;

    PalindromeDeletionModel &operator=(const PalindromeDeletionModel &) = default;

    PalindromeDeletionModel &operator=(PalindromeDeletionModel &&) = default;

    virtual ~PalindromeDeletionModel() = default;

    const DeletionTableVector &GetDeletionTable() const { return deletion_table_; }

    void SetDeletionTable(const DeletionTableVector &deletion_table) {
        deletion_table_ = deletion_table;
    }

    const germline_utils::ImmuneGeneDatabase* GetIgGeneDatabasePtr() const {
        return ig_gene_database_;
    }

    const germline_utils::ImmuneGeneDatabase &GetIgGeneDatabase() const {
        assert(ig_gene_database_ != nullptr);
        return *ig_gene_database_;
    }

    const std::vector<int> &GetDeletionLength() const { return deletion_length_; }

    using iterator = DeletionTableVector::iterator;
    using citerator = DeletionTableVector::const_iterator;

    iterator  begin ()       { return deletion_table_.begin(); }
    citerator begin () const { return deletion_table_.begin(); }
    citerator cbegin() const { return deletion_table_.cbegin(); }
    iterator  end   ()       { return deletion_table_.end(); }
    citerator end   () const { return deletion_table_.end(); }
    citerator cend  () const { return deletion_table_.cend(); }

    size_t size() const { return deletion_table_.size(); }

    PalindromeDeletionModel(std::ifstream &, const germline_utils::ImmuneGeneDatabase*);

    PalindromeDeletionModel(std::ifstream &, const germline_utils::ImmuneGeneDatabase &);

    double GetDeletionProbability(const size_t &gene_id, const int deletion_length) const;
};

std::ostream &operator<<(std::ostream &, const PalindromeDeletionModel &);

/***************************************************************************************************/

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

    const germline_utils::ChainDatabase* HC_database;

public:
    HCProbabilityRecombinationModel() = delete;
    HCProbabilityRecombinationModel(const HCProbabilityRecombinationModel &) = default;
    HCProbabilityRecombinationModel(HCProbabilityRecombinationModel &&) = default;
    HCProbabilityRecombinationModel &operator=(const HCProbabilityRecombinationModel &) = default;
    HCProbabilityRecombinationModel &operator=(HCProbabilityRecombinationModel &&) = default;

    virtual ~HCProbabilityRecombinationModel() = default;

    const IgGeneProbabilityModel &GetVGeneProbabilityModel() const {
        return V_gene_probability_model_;
    }

    const IgGeneProbabilityModel &GetDGeneProbabilityModel() const {
        return D_gene_probability_model_;
    }

    const IgGeneProbabilityModel &GetJGeneProbabilityModel() const {
        return J_gene_probability_model_;
    }

    const NongenomicInsertionModel &GetVDNongenomicInsertionModel() const {
        return VD_nongenomic_insertion_model_;
    }

    const NongenomicInsertionModel &GetDJNongenomicInsertionModel() const {
        return DJ_nongenomic_insertion_model_;
    }

    const PalindromeDeletionModel &GetVPalindromeDeletionModel() const {
        return V_palindrome_deletion_model_;
    }

    const PalindromeDeletionModel &GetJPalindromeDeletionModel() const {
        return J_palindrome_deletion_model_;
    }

    const PalindromeDeletionModel &GetDLeftPalindromeDeletionModel() const {
        return DLeft_palindrome_deletion_model_;
    }

    const PalindromeDeletionModel &GetDRightPalindromeDeletionModel() const {
        return DRight_palindrome_deletion_model_;
    }

    void SetVGeneProbabilityModel(const IgGeneProbabilityModel V_gene_probability_model) {
        V_gene_probability_model_ = V_gene_probability_model;
    }

    void SetDGeneProbabilityModel(const IgGeneProbabilityModel D_gene_probability_model) {
        D_gene_probability_model_ = D_gene_probability_model;
    }

    void SetJGeneProbabilityModel(const IgGeneProbabilityModel J_gene_probability_model) {
        J_gene_probability_model_ = J_gene_probability_model;
    }

    HCProbabilityRecombinationModel(std::ifstream &, const germline_utils::ChainDatabase*);
    HCProbabilityRecombinationModel(std::ifstream &, const germline_utils::ChainDatabase &);

    //double GetProbabilityByGenId(const IgGene& ig_gene) const;
    double GetProbabilityByGenId(germline_utils::SegmentType ig_gene_type,
                                 const size_t &id) const;
};

std::ostream &operator<<(std::ostream &, const HCProbabilityRecombinationModel &);

} // End namespace vdj_labeler

#include "recombination_model.hpp"

#include <cassert>

#include <boost/tokenizer.hpp>

namespace vdj_labeler {

using boost::tokenizer;
using boost::escaped_list_separator;
using Tokenizer = tokenizer<escaped_list_separator<char>>;

IgGeneProbabilityModel::IgGeneProbabilityModel(const IgGeneProbabilityVector &ig_gene_probabilities,
                                               const IgGeneDatabasePtrConst &ig_gene_database) :
    ig_gene_probabilities_(ig_gene_probabilities),
    ig_gene_database_(ig_gene_database) { }

size_t IgGeneProbabilityModel::size() const { return ig_gene_probabilities_.size(); }

IgGeneProbabilityModel::IgGeneProbabilityModel(std::ifstream &in,
                                               const IgGeneDatabasePtrConst &ig_gene_database) :
    ig_gene_database_(ig_gene_database) {
    assert(in.is_open());
    std::string line;
    std::vector <std::string> parsed_vector;
    ig_gene_probabilities_.resize(ig_gene_database_->size());
    while (getline(in, line)) {
        if (line.empty())
            break;
        Tokenizer tokenizer(line);
        parsed_vector.assign(tokenizer.begin(), tokenizer.end());
        assert(parsed_vector.size() == 2);
        size_t index_of_current_gene = ig_gene_database->GetIndexByName(parsed_vector.front());
        ig_gene_probabilities_[index_of_current_gene] = std::stod(parsed_vector.back());
    }
}

IgGeneProbabilityModel::IgGeneProbabilityModel(std::ifstream &in, const germline_utils::ImmuneGeneDatabase &database) :
    IgGeneProbabilityModel(in, std::make_shared<const germline_utils::ImmuneGeneDatabase>(database)) { }

std::ostream &operator<<(std::ostream &out, const IgGeneProbabilityModel &model) {
    for (size_t i = 0; i < model.size(); ++i)
        out << "Gene_id: " << model.GetIgGeneDatabase()->operator[](i).name()
            << ", " << "Gene probability: " << model.GetIgGeneProbabilities().at(i) << "\n";
    return out;
}

double IgGeneProbabilityModel::GetProbabilityByGenId(const size_t &id) const {
    assert(ig_gene_probabilities_.size() > id);
    return ig_gene_probabilities_[id];
}

/************************************************************************************************/

NongenomicInsertionModel::NongenomicInsertionModel(
    std::vector<double> insertion_probabilities,
    NongenomicInsertionMatrix transition_matrix) :
    insertion_probabilities_(insertion_probabilities),
    transition_matrix_(transition_matrix) { }

NongenomicInsertionModel::NongenomicInsertionModel(std::ifstream &in) {
    assert(in.is_open());
    std::string line;
    std::vector <std::string> parsed_vector;
    while (getline(in, line)) {
        if (line.empty())
            break;
        Tokenizer tokenizer(line);
        parsed_vector.assign(tokenizer.begin(), tokenizer.end());
        assert(parsed_vector.size() == 2);
        insertion_probabilities_.push_back(std::stod(parsed_vector[1]));
    }

    while (getline(in, line)) {
        if (line.empty())
            break;
        Tokenizer tokenizer(line);
        transition_matrix_.emplace_back(std::distance(tokenizer.begin(), tokenizer.end()));
        std::transform(tokenizer.begin(), tokenizer.end(),
                       transition_matrix_.back().begin(),
                       [](std::string str) { return std::stod(str); });
    }
}

double NongenomicInsertionModel::GetInsertionProbabilityByLength(
    const long unsigned int ins_length) const {
    assert(ins_length < insertion_probabilities_.size());
    return insertion_probabilities_[ins_length];
}

double NongenomicInsertionModel::GetTransitionProbability(char in, char out) const {
    using Alphabet = seqan::Dna5;
    size_t in_pos = static_cast<size_t>(seqan::ordValue(Alphabet(in)));
    size_t out_pos = static_cast<size_t>(seqan::ordValue(Alphabet(out)));
    assert(in_pos < 4);
    assert(out_pos < 4);
    return transition_matrix_[in_pos][out_pos];
}


std::ostream &operator<<(std::ostream &out, const NongenomicInsertionModel &model) {
    out << "Length, Probability" << "\n";
    for (auto it = model.GetInsertionProbabilities().cbegin();
         it != model.GetInsertionProbabilities().cend(); ++it)
        out << it - model.GetInsertionProbabilities().cbegin() << ", " << *it << "\n";
    out << "\n";

    out << "Markov chain transition matrix" << "\n";
    for (auto it1 = model.GetTransitionMatrix().cbegin();
         it1 != model.GetTransitionMatrix().cend(); ++it1) {
        for (auto it2 = it1->cbegin(); it2 != it1->cend(); ++it2)
            out << *it2 << " ";
        out << "\n";
    }
    return out;
}

double NongenomicInsertionModel::GetTransitionProbability(
    const std::pair<char, char> &transition) const {
    return GetTransitionProbability(transition.first, transition.second);
}

/**************************************************************************************************/

PalindromeDeletionModel::PalindromeDeletionModel(const DeletionTableVector &deletion_table,
                                                 const std::vector<int> &deletion_length,
                                                 const IgGeneDatabasePtrConst &ig_gene_database) :
    deletion_table_(deletion_table),
    deletion_length_(deletion_length),
    ig_gene_database_(ig_gene_database) { }


PalindromeDeletionModel::PalindromeDeletionModel(std::ifstream &in,
                                                 const IgGeneDatabasePtrConst &ig_gene_database) :
    ig_gene_database_(ig_gene_database) {
    assert(in.is_open());
    deletion_table_.resize(ig_gene_database->size());

    std::string line;
    std::vector <std::string> parsed_vector;
    getline(in, line);
    Tokenizer tokenizer(line);
    size_t palidrome_len_diversity = std::distance(tokenizer.begin(), tokenizer.end()) - 1;
    deletion_length_.resize(palidrome_len_diversity);

    auto tokenizer_it = tokenizer.begin();
    tokenizer_it++;

    for (auto it = deletion_length_.begin();
         it != deletion_length_.end();
         ++it, ++tokenizer_it) {
        *it = std::stoi(*tokenizer_it);
    }

    for (auto it = deletion_table_.begin(); it != deletion_table_.end(); ++it)
        it->resize(palidrome_len_diversity);

    while (getline(in, line)) {
        if (line.empty())
            break;
        Tokenizer tokenizer(line);
        parsed_vector.assign(tokenizer.begin(), tokenizer.end());
        assert(parsed_vector.size() == palidrome_len_diversity + 1);
        auto palidrome_len_it = 0;
        size_t index_of_current_gene = ig_gene_database->GetIndexByName(parsed_vector.front());
        for (auto parsed_it = parsed_vector.begin() + 1;
             parsed_it != parsed_vector.end();
             ++parsed_it, ++palidrome_len_it) {
            deletion_table_[index_of_current_gene][palidrome_len_it] = std::stod(*parsed_it);
        }
    }
}

PalindromeDeletionModel::PalindromeDeletionModel(std::ifstream &in,
                                                 const germline_utils::ImmuneGeneDatabase &ig_gene_database) :
    PalindromeDeletionModel(in, std::make_shared<const germline_utils::ImmuneGeneDatabase>(ig_gene_database)) { }

std::ostream &operator<<(std::ostream &out, const PalindromeDeletionModel &model) {
    if (model.GetDeletionTable().empty())
        return out;
    out << "Gene id, Length of palindrome: ";

    for (int length : model.GetDeletionLength())
        out << length << " ";
    out << "\n";

    for (size_t i = 0; i < model.size(); ++i) {
        out << model.GetIgGeneDatabase()->operator[](i).name() << " ";
        for (auto it = model.GetDeletionTable().at(i).begin();
             it != model.GetDeletionTable().at(i).end();
             ++it)
            out << *it << " ";
        out << "\n";
    }
    return out;
}

double PalindromeDeletionModel::GetDeletionProbability(const size_t &gene_id, const int deletion_length) const {
    auto it_deletion_length = std::find(deletion_length_.begin(), deletion_length_.end(), deletion_length);
    assert(it_deletion_length != deletion_length_.end());
    return deletion_table_[gene_id][it_deletion_length - deletion_length_.begin()];
}

/**************************************************************************************************/

HCProbabilityRecombinationModel::HCProbabilityRecombinationModel(std::ifstream &in,
                                                                 const HC_GenesDatabase_PtrConst &HC_db) :
    V_gene_probability_model_(in, std::make_shared<const germline_utils::ImmuneGeneDatabase>(HC_db->GetDb(germline_utils::VariableSegment))),
    D_gene_probability_model_(in, std::make_shared<const germline_utils::ImmuneGeneDatabase>(HC_db->GetDb(germline_utils::DiversitySegment))),
    J_gene_probability_model_(in, std::make_shared<const germline_utils::ImmuneGeneDatabase>(HC_db->GetDb(germline_utils::JoinSegment))),
    VD_nongenomic_insertion_model_(in),
    DJ_nongenomic_insertion_model_(in),
    V_palindrome_deletion_model_(in, std::make_shared<const germline_utils::ImmuneGeneDatabase>(HC_db->GetDb(germline_utils::VariableSegment))),
    J_palindrome_deletion_model_(in, std::make_shared<const germline_utils::ImmuneGeneDatabase>(HC_db->GetDb(germline_utils::JoinSegment))),
    DLeft_palindrome_deletion_model_(in, std::make_shared<const germline_utils::ImmuneGeneDatabase>(HC_db->GetDb(germline_utils::DiversitySegment))),
    DRight_palindrome_deletion_model_(in, std::make_shared<const germline_utils::ImmuneGeneDatabase>(HC_db->GetDb(germline_utils::DiversitySegment))),
    HC_database(HC_db)
{ }

HCProbabilityRecombinationModel::HCProbabilityRecombinationModel(std::ifstream &in,
                                                                 const germline_utils::ChainDatabase &HC_db) :
    HCProbabilityRecombinationModel(in, std::make_shared<germline_utils::ChainDatabase>(HC_db)) { }

std::ostream &operator<<(std::ostream &out, const HCProbabilityRecombinationModel &model) {
    out << "V gene probabilities:\n";
    out << model.GetVGeneProbabilityModel();
    out << "\nD gene probabilities:\n";
    out << model.GetDGeneProbabilityModel();
    out << "\nJ gene probabilities:\n";
    out << model.GetJGeneProbabilityModel();

   out << "\nVD nongenomic insertions model:\n";
   out << model.GetVDNongenomicInsertionModel();
   out << "\nDJ nongenomic insertions model:\n";
   out << model.GetDJNongenomicInsertionModel();

   out << "\nV palindrome deletion model:\n";
   out << model.GetVPalindromeDeletionModel();
   out << "\nJ palindrome deletion model:\n";
   out << model.GetJPalindromeDeletionModel();
   out << "\nD left palindrome deletion model:\n";
   out << model.GetDLeftPalindromeDeletionModel();
   out << "\nD right palindrome deletion model:\n";
   out << model.GetDRightPalindromeDeletionModel();
   return out;
}

double HCProbabilityRecombinationModel::GetProbabilityByGenId(germline_utils::SegmentType ig_gene_type,
                                                              const size_t &id) const {
    assert(ig_gene_type != germline_utils::SegmentType::VariableSegment);
    if (ig_gene_type == germline_utils::SegmentType::VariableSegment)
        return V_gene_probability_model_.GetProbabilityByGenId(id);
    if (ig_gene_type == germline_utils::SegmentType::DiversitySegment)
        return D_gene_probability_model_.GetProbabilityByGenId(id);
    // ig_gene_type == join_gene
    return J_gene_probability_model_.GetProbabilityByGenId(id);
}

}

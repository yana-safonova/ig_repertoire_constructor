//
// Created by Andrew Bzikadze on 3/9/17.
//

#include "shm_model.hpp"

#include <string>

#include <cassert>
#include <boost/tokenizer.hpp>

#include "verify.hpp"

using namespace shm_kmer_matrix_estimator;

namespace antevolo {

using boost::tokenizer;
using boost::escaped_list_separator;
using Tokenizer = tokenizer<escaped_list_separator<char>>;

ShmModel::ShmModel(const std::string &filename) {
    std::fstream in(filename);
    VERIFY_MSG(in.is_open(), std::string("File is not opened! Filename:") + filename);
    std::string line;
    std::vector<std::string> parsed_vector;

    getline(in, line); // skip the header

    auto string_to_bool = [](std::string &s) { return std::stoi(s) != 0; };

    while (getline(in, line)) {
        if (line.empty()) {
            break;
        }
        Tokenizer tokenizer(line);
        parsed_vector.assign(tokenizer.begin(), tokenizer.end());
        VERIFY_MSG(parsed_vector.size() == 23,
                   std::string("Length of parsed_vector = ") + std::to_string(parsed_vector.size()) +
                   ".Correct = 23");

        size_t i = 1;
        beta_fr_params_.push_back({std::stod(parsed_vector[i]),
                                   std::stod(parsed_vector[i + 1])});
        i += 2;
        beta_cdr_params_.push_back({std::stod(parsed_vector[i]),
                                    std::stod(parsed_vector[i + 1])});
        i += 2;
        beta_full_params_.push_back({std::stod(parsed_vector[i]),
                                     std::stod(parsed_vector[i + 1])});
        i += 2;

        dirichlet_params_.push_back({std::stod(parsed_vector[i]),
                                     std::stod(parsed_vector[i + 1]),
                                     std::stod(parsed_vector[i + 2])});
        i += 3;

        beta_fr_success_mle_.push_back(  string_to_bool(parsed_vector[i++]));
        beta_cdr_success_mle_.push_back( string_to_bool(parsed_vector[i++]));
        beta_full_success_mle_.push_back(string_to_bool(parsed_vector[i++]));

        dirichlet_success_mle_.push_back(string_to_bool(parsed_vector[i++]));


        start_point_beta_fr_params_.push_back({std::stod(parsed_vector[i]),
                                               std::stod(parsed_vector[i + 1])});
        i += 2;
        start_point_beta_cdr_params_.push_back({std::stod(parsed_vector[i]),
                                                std::stod(parsed_vector[i + 1])});
        i += 2;
        start_point_beta_full_params_.push_back({std::stod(parsed_vector[i]),
                                                 std::stod(parsed_vector[i + 1])});
        i += 2;
        start_point_dirichlet_params_.push_back({std::stod(parsed_vector[i]),
                                                 std::stod(parsed_vector[i + 1]),
                                                 std::stod(parsed_vector[i + 2])});
        i += 3;
    }
}

size_t ShmModel::size() const {
    return beta_fr_params_.size();
}

std::ostream &operator<<(std::ostream &os, const ShmModel &obj) {
    size_t kmer_len = obj.beta_fr_params().kmer_len();
    std::vector<std::string> kmers(KmerUtils::GenerateAllKmersFixedLength(kmer_len));
    const char sep = ',';
    os << sep <<
       "beta_FR_shape1"               << sep << "beta_FR_shape2"               << sep <<
       "beta_CDR_shape1"              << sep << "beta_CDR_shape2"              << sep <<
       "beta_FULL_shape1"             << sep << "beta_FULL_shape2"             << sep <<
       "dir_shape1"                   << sep << "dir_shape2"                   << sep << "dir_shape3" << sep <<

       "success_optim_beta_FR"        << sep <<
       "success_optim_beta_CDR"       << sep <<
       "success_optim_beta_FULL"      << sep <<
       "success_optim_dir"            << sep <<

       "start_point_beta_FR_shape1"   << sep << "start_point_beta_FR_shape2"   << sep <<
       "start_point_beta_CDR_shape1"  << sep << "start_point_beta_CDR_shape2"  << sep <<
       "start_point_beta_FULL_shape1" << sep << "start_point_beta_FULL_shape2" << sep <<
       "start_point_dir_shape1"       << sep << "start_point_dir_shape2"       << sep << "start_point_dir_shape3" <<
       "\n";
    for (size_t i = 0; i < obj.size(); ++i) {
        os << kmers[i] << sep <<
            obj.beta_fr_params()[i][0]               << sep << obj.beta_fr_params()[i][1]   << sep <<
            obj.beta_cdr_params()[i][0]              << sep << obj.beta_cdr_params()[i][1]  << sep <<
            obj.beta_full_params()[i][0]             << sep << obj.beta_full_params()[i][1] << sep <<
            obj.dirichlet_params()[i][0]             << sep << obj.dirichlet_params()[i][1] << sep <<
                obj.dirichlet_params()[i][2]         << sep <<

            obj.beta_fr_success_mle()[i]             << sep <<
            obj.beta_cdr_success_mle()[i]            << sep <<
            obj.beta_full_success_mle()[i]           << sep <<
            obj.dirichlet_success_mle()[i]           << sep <<

            obj.start_point_beta_fr_params()[i][0]   << sep << obj.start_point_beta_fr_params()[i][1]   << sep <<
            obj.start_point_beta_cdr_params()[i][0]  << sep << obj.start_point_beta_cdr_params()[i][1]  << sep <<
            obj.start_point_beta_full_params()[i][0] << sep << obj.start_point_beta_full_params()[i][1] << sep <<
            obj.start_point_dirichlet_params()[i][0] << sep << obj.start_point_dirichlet_params()[i][1] << sep <<
                obj.start_point_dirichlet_params()[i][2] <<
            '\n';

    }
    return os;
}

} // End namespace antevolo

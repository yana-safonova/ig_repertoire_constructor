//
// Created by Andrew Bzikadze on 3/9/17.
//

#include "shm_model.hpp"

#include <string>

#include <cassert>
#include <boost/tokenizer.hpp>

#include "verify.hpp"

using namespace antevolo;

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

} // End namespace antevolo

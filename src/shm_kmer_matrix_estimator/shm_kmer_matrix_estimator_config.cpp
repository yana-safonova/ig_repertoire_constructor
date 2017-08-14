//
// Created by Andrew Bzikadze on 5/14/16.
//

#include <stdexcept>
#include <algorithm>

#include "shm_kmer_matrix_estimator_config.hpp"

#include "../include/config_common.hpp"

namespace shm_kmer_matrix_estimator {

std::string error_message_strategy(const std::string &what_about,
                                   const std::string &supplied_method,
                                   const std::vector<std::string> &available_methods,
                                   const std::string &where = "shm_kmer_matrix_estimator_config") {
    std::string message(":::(");
    message += where + ") wrong ";
    message += what_about;
    message += " method: ";
    message += supplied_method;
    message += ". Available: ";
    for (auto available_method = available_methods.begin();
         available_method != available_methods.end();
         ++available_method) {
        message += *available_method + ", ";
    }
    message += ":::";
    return message;
}

// IO parameters START
void updateIO(shm_kmer_matrix_estimator_config::io_params &io) {
    io.output.log_filename = path::append_path(io.output.output_dir, io.output.log_filename);
    io.output.output_filename_fr = path::append_path(io.output.output_dir, io.output.output_filename_fr);
    io.output.output_filename_cdr = path::append_path(io.output.output_dir, io.output.output_filename_cdr);
}

void load(shm_kmer_matrix_estimator_config::io_params::input_params &input_params, boost::property_tree::ptree const &pt, bool) {
    using config_common::load;
    load(input_params.v_alignments, pt, "v_alignments");
    load(input_params.cdr_details, pt, "cdr_details");
}

void load(shm_kmer_matrix_estimator_config::io_params::output_params &output, boost::property_tree::ptree const &pt, bool) {
    using config_common::load;
    load(output.output_dir, pt, "output_dir");
    load(output.log_filename, pt, "log_filename");
    load(output.output_filename_fr, pt, "output_filename_fr");
    load(output.output_filename_cdr, pt, "output_filename_cdr");
}

void load(shm_kmer_matrix_estimator_config::io_params &io, boost::property_tree::ptree const &pt, bool) {
    using config_common::load;
    load(io.input, pt, "input_params");
    load(io.output, pt, "output_params");
    updateIO(io);
}
// IO parameters FINISH

// Alignment Checker parameters START
void load(shm_kmer_matrix_estimator_config::alignment_checker_params &achp,
          boost::property_tree::ptree const &pt, bool) {
    using config_common::load;
    using AlignmentCheckerMethod = shm_kmer_matrix_estimator_config::alignment_checker_params::AlignmentCheckerMethod;
    std::string method_str(pt.get<std::string>("alignment_checker_method"));
    std::string method_str_lowercase(method_str);
    std::transform(method_str.begin(), method_str.end(),
                   method_str_lowercase.begin(), ::tolower);
    if (method_str_lowercase == "nogaps") {
        achp.alignment_checker_method = AlignmentCheckerMethod::NoGaps;
    } else {
        std::string message = error_message_strategy("alignment checker strategy",
                                                     method_str,
                                                     achp.alignment_checker_method_names);
        throw std::invalid_argument(message);
    }

    using FunctionalityMethod = shm_kmer_matrix_estimator_config::alignment_checker_params::FunctionalityMethod;
    std::string functionality_method(pt.get<std::string>("functionality_method"));
    if (functionality_method == "all") {
        achp.functionality_method = FunctionalityMethod::all;
    } else if (functionality_method == "productive") {
        achp.functionality_method = FunctionalityMethod::productive;
    } else if (functionality_method == "nonproductive") {
        achp.functionality_method = FunctionalityMethod::nonproductive;
    } else {
        std::string message = error_message_strategy("alignment checker functionality checker",
                                                     functionality_method,
                                                     achp.functionality_methods_names);
        throw std::invalid_argument(message);
    }
}
// Alignment Checker parameters FINISH

// Alignment Cropper parameters START
void load(shm_kmer_matrix_estimator_config::alignment_cropper_params::upto_reliable_kmer_cropper_params &rkmp,
          boost::property_tree::ptree const &pt, bool) {
    using config_common::load;
    load(rkmp.kmer_len, pt, "kmer_len");
    load(rkmp.hash_base, pt, "hash_base");
}

void load(shm_kmer_matrix_estimator_config::alignment_cropper_params &acrp,
          boost::property_tree::ptree const &pt, bool) {
    using config_common::load;
    using AlignmentCropperMethod = shm_kmer_matrix_estimator_config::alignment_cropper_params::AlignmentCropperMethod;
    std::string method_str(pt.get<std::string>("alignment_cropper_method"));
    std::string method_str_lowercase(method_str);
    std::transform(method_str.begin(), method_str.end(),
                   method_str_lowercase.begin(), ::tolower);
    if (method_str_lowercase == "uptolastreliablekmer") {
        acrp.alignment_cropper_method = AlignmentCropperMethod::UptoLastReliableKMer;
        load(acrp.rkmp, pt, "method_params");
    } else {
        std::string message = error_message_strategy("alignment cropper",
                                                     method_str,
                                                     acrp.alignment_cropper_method_names);
        throw std::invalid_argument(message);
    }
}
// Alignment Cropper parameters FINISH

// Mismatch Finder parameters START
void load(shm_kmer_matrix_estimator_config::mutations_strategy_params::trivial_mutations_strategy_params &,
          boost::property_tree::ptree const &, bool) {}

void load(shm_kmer_matrix_estimator_config::mutations_strategy_params::no_kneighbours_mutations_strategy_params &,
          boost::property_tree::ptree const &, bool) {}

void load(shm_kmer_matrix_estimator_config::mutations_strategy_params &mfp,
          boost::property_tree::ptree const &pt, bool) {
    using config_common::load;
    using MutationsStrategyMethod = shm_kmer_matrix_estimator_config::mutations_strategy_params::MutationsStrategyMethod;
    std::string method_str(pt.get<std::string>("mutations_strategy_method"));
    std::string method_str_lowercase(method_str);
    std::transform(method_str.begin(), method_str.end(),
                   method_str_lowercase.begin(), ::tolower);
    if (method_str_lowercase == "trivial") {
        mfp.mutations_strategy_method = MutationsStrategyMethod::Trivial;
        load(mfp.tmfp, pt, "method_params");
    } else if (method_str_lowercase == "nokneighbours") {
        mfp.mutations_strategy_method = MutationsStrategyMethod::NoKNeighbours;
        load(mfp.nknmfp, pt, "method_params");
    } else {
        std::string message = error_message_strategy("mutation strategy",
                                                     method_str,
                                                     mfp.mutation_strategy_method_names);
        throw std::invalid_argument(message);
    }
    load(mfp.kmer_len, pt, "kmer_len");
}
// Mismatch Finder parameters FINISH

// Main load functions
void load(shm_kmer_matrix_estimator_config &cfg, boost::property_tree::ptree const &pt, bool complete) {
    using config_common::load;
    load(cfg.io, pt, "io_params", complete);
    load(cfg.achp, pt, "alignment_checker_params");
    load(cfg.acrp, pt, "alignment_cropper_params");
    load(cfg.mfp, pt, "mutations_strategy_params");
}

void load(shm_kmer_matrix_estimator_config &cfg, std::string const &filename) {
    boost::property_tree::ptree pt;
    boost::property_tree::read_info(filename, pt);
    load(cfg, pt, true);
}

std::istream &operator>>(std::istream &in, shm_kmer_matrix_estimator_config::mutations_strategy_params::MutationsStrategyMethod &strategy) {
    std::string token_original;
    in >> token_original;
    std::string token(token_original.size(), ' ');
    std::transform(token_original.begin(), token_original.end(),
                   token.begin(), ::tolower);
    using MutationStrategyMethod = shm_kmer_matrix_estimator_config::mutations_strategy_params::MutationsStrategyMethod;
    if (token == "trivial")
        strategy = MutationStrategyMethod::Trivial;
    else if (token == "nokneighbours")
        strategy = MutationStrategyMethod::NoKNeighbours;
    else {
        std::string message = error_message_strategy("mutation strategy",
                                                     token_original,
                                                     shm_kmer_matrix_estimator_config::mutations_strategy_params::
                                                     mutation_strategy_method_names,
                                                     "command line argument");
        throw std::invalid_argument(message);
    }
    return in;
}

const std::vector<std::string> shm_kmer_matrix_estimator_config::alignment_checker_params::alignment_checker_method_names =
    {std::string("NoGaps")};
const std::vector<std::string> shm_kmer_matrix_estimator_config::alignment_checker_params::functionality_methods_names =
    {std::string("all"), std::string("productive"), std::string("nonproductive")};
const std::vector<std::string> shm_kmer_matrix_estimator_config::alignment_cropper_params::alignment_cropper_method_names =
    {std::string("UptoLastReliableKMer")};
const std::vector<std::string> shm_kmer_matrix_estimator_config::mutations_strategy_params::mutation_strategy_method_names =
    {std::string("Trivial"), std::string("NoKNeighbours")};

} // End namespace shm_kmer_matrix_estimator

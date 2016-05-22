//
// Created by Andrew Bzikadze on 5/14/16.
//

#include "shm_config.hpp"

#include "../include/config_common.hpp"

// IO parameters START
void updateIO(shm_config::io_params &io) {
    io.output.log_filename = path::append_path(io.output.output_dir, io.output.log_filename);
    io.output.output_filename = path::append_path(io.output.output_dir, io.output.output_filename);
}

void load(shm_config::io_params::input_params &input_params, boost::property_tree::ptree const &pt, bool) {
    using config_common::load;
    load(input_params.input_filename, pt, "input_filename");
}

void load(shm_config::io_params::output_params &output, boost::property_tree::ptree const &pt, bool) {
    using config_common::load;
    load(output.output_dir, pt, "output_dir");
    load(output.log_filename, pt, "log_filename");
    load(output.output_filename, pt, "output_filename");
}

void load(shm_config::io_params& io, boost::property_tree::ptree const& pt, bool) {
    using config_common::load;
    load(io.input, pt, "input_params");
    load(io.output, pt, "output_params");
    updateIO(io);
}
// IO parameters FINISH

// Alignment Checker parameters START
void load(shm_config::alignment_checker_params& achp,
          boost::property_tree::ptree const& pt, bool) {
    using config_common::load;
    using AlignmentCheckerMethod = shm_config::alignment_checker_params::AlignmentCheckerMethod;
    if (pt.get<std::string>("alignment_checker_method") == "NoGapsAlignmentChecker") {
        achp.alignment_checker_method = AlignmentCheckerMethod::NoGapsAlignmentChecker;
    }
}
// Alignment Checker parameters FINISH

// Alignment Cropper parameters START
void load(shm_config::alignment_cropper_params::upto_reliable_kmer_cropper_params& rkmp,
          boost::property_tree::ptree const& pt, bool) {
    using config_common::load;
    load(rkmp.kmer_len, pt, "kmer_len");
    load(rkmp.hash_base, pt, "hash_base");
}

void load(shm_config::alignment_cropper_params& acrp,
          boost::property_tree::ptree const& pt, bool) {
    using config_common::load;
    using AlignmentCropperMethod = shm_config::alignment_cropper_params::AlignmentCropperMethod;
    if (pt.get<std::string>("alignment_cropper_method") == "UptoLastReliableKMer") {
        acrp.alignment_cropper_method = AlignmentCropperMethod::UptoLastReliableKMer;
        load(acrp.rkmp, pt, "method_params");
    }
}
// Alignment Cropper parameters FINISH

// Mismatch Finder parameters START
void load(shm_config::mutations_strategy_params::trivial_mutations_strategy_params& tmfp,
          boost::property_tree::ptree const&, bool) { }

void load(shm_config::mutations_strategy_params::no_kneighbours_mutations_strategy_params& nknmfp,
          boost::property_tree::ptree const& pt, bool) {
    using config_common::load;
    load(nknmfp.kmer_len, pt, "kmer_len");
}

void load(shm_config::mutations_strategy_params& mfp,
          boost::property_tree::ptree const& pt, bool) {
    using config_common::load;
    using MutationsStrategyMethod = shm_config::mutations_strategy_params::MutationsStrategyMethod;
    if (pt.get<std::string>("mutations_strategy_method") == "Trivial") {
        mfp.mutations_strategy_method = MutationsStrategyMethod::Trivial;
        load(mfp.tmfp, pt, "method_params");
    } else if (pt.get<std::string>("mutations_strategy_method") == "NoKNeighbours") {
        mfp.mutations_strategy_method = MutationsStrategyMethod::NoKNeighbours;
        load(mfp.nknmfp, pt, "method_params");
    }
}
// Mismatch Finder parameters FINISH

// Main load functions
void load(shm_config& cfg, boost::property_tree::ptree const& pt, bool complete) {
    using config_common::load;
    load(cfg.io, pt, "io_params", complete);
    load(cfg.achp, pt, "alignment_checker_params");
    load(cfg.acrp, pt, "alignment_cropper_params");
    load(cfg.mfp, pt, "mutations_strategy_params");
}

void load(shm_config& cfg, std::string const& filename) {
    boost::property_tree::ptree pt;
    boost::property_tree::read_info(filename, pt);
    load(cfg, pt, true);
}

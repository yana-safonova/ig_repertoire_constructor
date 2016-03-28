#include <logger/logger.hpp>
#include "vj_finder_config.hpp"
#include <config_common.hpp>
#include <boost/program_options.hpp>
#include <build_info.hpp>

namespace vj_finder {

    void load(vjf_config::run_params &rp, boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(rp.num_threads, pt, "num_threads");
    }

    void load(vjf_config::io_params::input_params& ip, boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(ip.input_reads, pt, "input_reads");
    }

    void load(vjf_config::io_params::output_params::output_details& od,
              boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(od.compress, pt, "compress");
        load(od.fix_spaces, pt, "fix_spaces");
        load(od.separator, pt, "separator");
        load(od.verbose, pt, "verbose");
    }

    void load(vjf_config::io_params::output_params::output_files& of,
              boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(of.add_info_filename, pt, "add_info_filename");
        load(of.bad_output_filename, pt, "bad_output_filename");
        load(of.discard_info_filename, pt, "discard_info_filename");
        load(of.output_dir, pt, "output_dir");
        load(of.output_filename, pt, "output_filename");
        load(of.valignments_filename, pt, "valignments_filename");
    }

    void load(vjf_config::io_params::output_params& op, boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(op.od, pt, "od");
        load(op.of, pt, "of");
    }

    void load(vjf_config::io_params &iop, boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(iop.input, pt, "input");
        load(iop.output, pt, "output");
    }

    void load(vjf_config::algorithm_params::aligner_params &ap, boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(ap.K, pt, "K");
        load(ap.left_uncovered_limit, pt, "left_uncovered_limit");
        load(ap.min_j_segment_length, pt, "min_j_segment_length");
        load(ap.min_k_coverage, pt, "min_k_coverage");
        load(ap.min_k_coverage_j, pt, "min_k_coverage_j");
        load(ap.min_v_segment_length, pt, "min_v_segment_length");
        load(ap.right_uncovered_limit, pt, "right_uncovered_limit");
        load(ap.word_size_j, pt, "word_size_j");
    }

    void load(vjf_config::algorithm_params::germline_params &gp, boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(gp.db_directory, pt, "db_directory");
        load(gp.loci, pt, "loci");
        load(gp.organism, pt, "organism");
        load(gp.pseudogenes, pt, "pseudogenes");
    }

    void load(vjf_config::algorithm_params::query_params &qp, boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(qp.max_candidates, pt, "max_candidates");
        load(qp.max_candidates_j, pt, "max_candidates_j");
        load(qp.fix_strand, pt, "fix_strand");
        load(qp.consistent_loci, pt, "consistent_loci");
        load(qp.min_len, pt, "min_len");
    }

    void load(vjf_config::algorithm_params::fix_crop_fill_params &fxp, boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(fxp.crop_left, pt, "crop_left");
        load(fxp.crop_right, pt, "crop_right");
        load(fxp.fill_left, pt, "fill_left");
        load(fxp.fill_right, pt, "fill_right");
        load(fxp.fix_left, pt, "fix_left");
        load(fxp.fix_right, pt, "fix_right");
    }

    void load(vjf_config::algorithm_params::alignment_scoring_params &asp, boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(asp.gap_extention_cost, pt, "gap_extention_cost");
        load(asp.gap_opening_cost, pt, "gap_opening_cost");
        load(asp.max_global_gap, pt, "max_global_gap");
        load(asp.max_local_deletions, pt, "max_local_deletions");
        load(asp.max_local_insertions, pt, "max_local_insertions");
    }

    void load(vjf_config::algorithm_params &algop, boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(algop.ap, pt, "ap");
        load(algop.gp, pt, "gp");
        load(algop.qp, pt, "qp");
        load(algop.fxp, pt, "fxp");
        load(algop.asp, pt, "asp");
    }

    void load(vjf_config &cfg, boost::property_tree::ptree const &pt, bool complete) {
        using config_common::load;
        load(cfg.rp, pt, "rp", complete);
        load(cfg.iop, pt, "iop", complete);
        load(cfg.algop, pt, "algop", complete);
    }

    void load(vjf_config &cfg, std::string const &filename) {
        boost::property_tree::ptree pt;
        boost::property_tree::read_info(filename, pt);
        load(cfg, pt, true);
    }
}
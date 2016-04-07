#include <logger/logger.hpp>
#include "vj_finder_config.hpp"
#include <config_common.hpp>
#include <boost/program_options.hpp>
#include <build_info.hpp>

namespace vj_finder {

    void load(vjf_config::RunParams &rp, boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(rp.num_threads, pt, "num_threads");
    }

    void load(vjf_config::IOParams::InputParams::GermlineInput &gi, boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(gi.germline_filenames_config, pt, "germline_filenames_config");
        load(gi.ig_dir, pt, "ig_dir");
        load(gi.tcr_dir, pt, "tcr_dir");
    }

    void update_input_config(vjf_config::IOParams::InputParams & ip) {
        ip.germline_input.germline_filenames_config = path::append_path(ip.config_dir,
                                                                        ip.germline_input.germline_filenames_config);
    }

    void load(vjf_config::IOParams::InputParams & ip, boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(ip.input_reads, pt, "input_reads");
        load(ip.config_dir, pt, "config_dir");
        load(ip.germline_input, pt, "germline_input");
        update_input_config(ip);
    }

    void load(vjf_config::IOParams::OutputParams::OutputDetails & od,
              boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(od.compress, pt, "compress");
        load(od.fix_spaces, pt, "fix_spaces");
        load(od.separator, pt, "separator");
        load(od.verbose, pt, "verbose");
    }

    void update_output_files_config(vjf_config::IOParams::OutputParams::OutputFiles & of) {
        of.add_info_filename = path::append_path(of.output_dir, of.add_info_filename);
        of.bad_output_filename = path::append_path(of.output_dir, of.bad_output_filename);
        of.discard_info_filename = path::append_path(of.output_dir, of.discard_info_filename);
        of.output_filename = path::append_path(of.output_dir, of.output_filename);
        of.valignments_filename = path::append_path(of.output_dir, of.valignments_filename);
    }

    void load(vjf_config::IOParams::OutputParams::OutputFiles & of,
              boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(of.add_info_filename, pt, "add_info_filename");
        load(of.bad_output_filename, pt, "bad_output_filename");
        load(of.discard_info_filename, pt, "discard_info_filename");
        load(of.output_dir, pt, "output_dir");
        load(of.output_filename, pt, "output_filename");
        load(of.valignments_filename, pt, "valignments_filename");
        update_output_files_config(of);
    }

    void load(vjf_config::IOParams::OutputParams & op, boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(op.output_details, pt, "output_details");
        load(op.output_files, pt, "output_files");
    }

    void load(vjf_config::IOParams &iop, boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(iop.input_params, pt, "input_params");
        load(iop.output_params, pt, "output_params");
    }

    void load(vjf_config::AlgorithmParams::AlignerParams &ap, boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(ap.word_size_v, pt, "word_size_v");
        load(ap.left_uncovered_limit, pt, "left_uncovered_limit");
        load(ap.min_j_segment_length, pt, "min_j_segment_length");
        load(ap.min_k_coverage_v, pt, "min_k_coverage_v");
        load(ap.min_k_coverage_j, pt, "min_k_coverage_j");
        load(ap.min_v_segment_length, pt, "min_v_segment_length");
        load(ap.right_uncovered_limit, pt, "right_uncovered_limit");
        load(ap.word_size_j, pt, "word_size_j");
    }

    void load(vjf_config::AlgorithmParams::GermlineParams &gp, boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(gp.germline_dir, pt, "germline_dir");
        load(gp.loci, pt, "loci");
        load(gp.organism, pt, "organism");
        load(gp.pseudogenes, pt, "pseudogenes");
    }

    void load(vjf_config::AlgorithmParams::QueryParams &qp, boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(qp.max_candidates_v, pt, "max_candidates_v");
        load(qp.max_candidates_j, pt, "max_candidates_j");
        load(qp.fix_strand, pt, "fix_strand");
        load(qp.min_len, pt, "min_len");
    }

    void load(vjf_config::AlgorithmParams::FixCropFillParams &fxp, boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(fxp.crop_left, pt, "crop_left");
        load(fxp.crop_right, pt, "crop_right");
        load(fxp.fill_left, pt, "fill_left");
        load(fxp.fill_right, pt, "fill_right");
        load(fxp.fix_left, pt, "fix_left");
        load(fxp.fix_right, pt, "fix_right");
    }

    void load(vjf_config::AlgorithmParams::ScoringParams::VScoringParams &vs,
              boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(vs.gap_extention_cost, pt, "gap_extention_cost");
        load(vs.gap_opening_cost, pt, "gap_opening_cost");
        load(vs.match_reward, pt, "match_reward");
        load(vs.max_global_gap, pt, "max_global_gap");
        load(vs.max_local_deletions, pt, "max_local_deletions");
        load(vs.max_local_insertions, pt, "max_local_insertions");
        load(vs.mismatch_extention_cost, pt, "mismatch_extention_cost");
        load(vs.mismatch_opening_cost, pt, "mismatch_opening_cost");
    }

    void load(vjf_config::AlgorithmParams::ScoringParams::JScoringParams &js,
              boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(js.gap_extention_cost, pt, "gap_extention_cost");
        load(js.gap_opening_cost, pt, "gap_opening_cost");
        load(js.match_reward, pt, "match_reward");
        load(js.max_global_gap, pt, "max_global_gap");
        load(js.max_local_deletions, pt, "max_local_deletions");
        load(js.max_local_insertions, pt, "max_local_insertions");
        load(js.mismatch_extention_cost, pt, "mismatch_extention_cost");
        load(js.mismatch_opening_cost, pt, "mismatch_opening_cost");
    }

    void load(vjf_config::AlgorithmParams::ScoringParams &asp, boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(asp.v_scoring, pt, "v_scoring");
        load(asp.j_scoring, pt, "j_scoring");
    }

    void load(vjf_config::AlgorithmParams &algop, boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(algop.aligner_params, pt, "aligner_params");
        load(algop.germline_params, pt, "germline_params");
        load(algop.query_params, pt, "query_params");
        load(algop.fix_crop_fill_params, pt, "fix_crop_fill_params");
        load(algop.scoring_params, pt, "scoring_params");
    }

    void load(vjf_config &cfg, boost::property_tree::ptree const &pt, bool complete) {
        using config_common::load;
        load(cfg.run_params, pt, "run_params", complete);
        load(cfg.io_params, pt, "io_params", complete);
        load(cfg.algorithm_params, pt, "algorithm_params", complete);
    }

    void load(vjf_config &cfg, std::string const &filename) {
        boost::property_tree::ptree pt;
        boost::property_tree::read_info(filename, pt);
        load(cfg, pt, true);
    }
}
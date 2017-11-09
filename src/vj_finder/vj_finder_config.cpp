#include <logger/logger.hpp>
#include "vj_finder_config.hpp"
#include <config_common.hpp>
#include <boost/program_options.hpp>
#include <build_info.hpp>
#include "germline_utils/germline_config.hpp"

namespace vj_finder {

    void load(VJFinderConfig::RunParams &rp, boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(rp.num_threads, pt, "num_threads");
    }

    void update_input_config(VJFinderConfig::IOParams::InputParams & ip) {
        ip.germline_input.germline_filenames_config = path::append_path(ip.config_dir,
                                                                        ip.germline_input.germline_filenames_config);
    }

    void load(VJFinderConfig::IOParams::InputParams & ip, boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(ip.input_reads, pt, "input_reads");
        load(ip.config_dir, pt, "config_dir");
        load(ip.germline_input, pt, "germline_input");
        //update_input_config(ip);
    }

    void load(VJFinderConfig::IOParams::OutputParams::OutputDetails & od, boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(od.fix_spaces, pt, "fix_spaces");
        load(od.num_aligned_candidates, pt, "num_aligned_candidates");
        std::string columns_str;
        load(od.alignment_columns, pt, "alignment_columns");
    }

    void update_output_files_config(VJFinderConfig::IOParams::OutputParams::OutputFiles & of) {
        of.log_filename = path::append_path(of.output_dir, of.log_filename);
        of.alignment_info_fname = path::append_path(of.output_dir, of.alignment_info_fname);
        of.filtered_reads_fname = path::append_path(of.output_dir, of.filtered_reads_fname);
        of.filtering_info_filename = path::append_path(of.output_dir, of.filtering_info_filename);
        of.cleaned_reads_fname = path::append_path(of.output_dir, of.cleaned_reads_fname);
        of.valignments_filename = path::append_path(of.output_dir, of.valignments_filename);
    }

    void load(VJFinderConfig::IOParams::OutputParams::OutputFiles & of, boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(of.log_filename, pt, "log_filename");
        load(of.alignment_info_fname, pt, "alignment_info_fname");
        load(of.filtered_reads_fname, pt, "filtered_reads_fname");
        load(of.filtering_info_filename, pt, "filtering_info_filename");
        load(of.output_dir, pt, "output_dir");
        load(of.cleaned_reads_fname, pt, "cleaned_reads_fname");
        load(of.valignments_filename, pt, "valignments_filename");
        //update_output_files_config(of);
    }

    void load(VJFinderConfig::IOParams::OutputParams & op, boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(op.output_details, pt, "output_details");
        load(op.output_files, pt, "output_files");
    }

    void load(VJFinderConfig::IOParams &iop, boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(iop.input_params, pt, "input_params");
        load(iop.output_params, pt, "output_params");
    }

    void load(VJFinderConfig::AlgorithmParams::AlignerParams &ap, boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(ap.word_size_v, pt, "word_size_v");
        load(ap.word_size_j, pt, "word_size_j");
        load(ap.min_k_coverage_v, pt, "min_k_coverage_v");
        load(ap.min_k_coverage_j, pt, "min_k_coverage_j");
        load(ap.max_candidates_v, pt, "max_candidates_v");
        load(ap.max_candidates_j, pt, "max_candidates_j");
        load(ap.fix_strand, pt, "fix_strand");
    }

    void load(VJFinderConfig::AlgorithmParams::FilteringParams &fp, boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(fp.enable_filtering, pt, "enable_filtering");
        load(fp.left_uncovered_limit, pt, "left_uncovered_limit");
        load(fp.right_uncovered_limit, pt, "right_uncovered_limit");
        load(fp.min_j_segment_length, pt, "min_j_segment_length");
        load(fp.min_v_segment_length, pt, "min_v_segment_length");
        load(fp.min_aligned_length, pt, "min_aligned_length");
    }

    VJFinderConfig::AlgorithmParams::FixCropFillParams::FixCropFillAlgorithm get_fcf_algorithm(std::string str) {
        if(str == "aggressive_fcf")
            return VJFinderConfig::AlgorithmParams::FixCropFillParams::FixCropFillAlgorithm::AggressiveFCFAlgorithm;
        VERIFY_MSG(false, "FCF algorithm was not recognized");
        return VJFinderConfig::AlgorithmParams::FixCropFillParams::FixCropFillAlgorithm::UnknowmFCFAlgorithm;
    }

    void load(VJFinderConfig::AlgorithmParams::FixCropFillParams &fxp, boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(fxp.enable_fixing, pt, "enable_fixing");
        load(fxp.enable_filling, pt, "enable_filling");
        load(fxp.enable_cropping, pt, "enable_cropping");
        load(fxp.fix_left, pt, "fix_left");
        load(fxp.fix_right, pt, "fix_right");
        load(fxp.crop_left, pt, "crop_left");
        load(fxp.crop_right, pt, "crop_right");
        load(fxp.fill_left, pt, "fill_left");
        load(fxp.fill_right, pt, "fill_right");
        std::string tmp;
        load(tmp, pt, "fcf_algorithm");
        fxp.fcf_algorithm = get_fcf_algorithm(tmp);
    }

    void load(VJFinderConfig::AlgorithmParams::ScoringParams::VScoringParams &vs,
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

    void load(VJFinderConfig::AlgorithmParams::ScoringParams::JScoringParams &js,
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

    void load(VJFinderConfig::AlgorithmParams::ScoringParams &asp, boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(asp.v_scoring, pt, "v_scoring");
        load(asp.j_scoring, pt, "j_scoring");
    }

    void load(VJFinderConfig::AlgorithmParams &algop, boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(algop.aligner_params, pt, "aligner_params");
        load(algop.germline_params, pt, "germline_params");
        load(algop.filtering_params, pt, "filtering_params");
        load(algop.fix_crop_fill_params, pt, "fix_crop_fill_params");
        load(algop.scoring_params, pt, "scoring_params");
    }

    void load(VJFinderConfig &cfg, boost::property_tree::ptree const &pt, bool complete) {
        using config_common::load;
        load(cfg.run_params, pt, "run_params", complete);
        load(cfg.io_params, pt, "io_params", complete);
        load(cfg.algorithm_params, pt, "algorithm_params", complete);
    }

    void load(VJFinderConfig &cfg, std::string const &filename) {
        boost::property_tree::ptree pt;
        boost::property_tree::read_info(filename, pt);
        load(cfg, pt, true);
    }

    using OutputDetails = VJFinderConfig::IOParams::OutputParams::OutputDetails;

    using VJFAlignmentInfoColumnTypeEnum = OutputDetails::AlignmentInfoColumnTypes::ColumnTypeEnum;

    const std::map<std::string, VJFAlignmentInfoColumnTypeEnum> OutputDetails::AlignmentInfoColumnTypes::string_to_column_type =
            std::map<std::string, VJFAlignmentInfoColumnTypeEnum> {
                    {"Read_name", ReadName},
                    {"Chain_type", ChainType},
                    {"V_hit", VHit},
                    {"V_start_pos", VStartPos},
                    {"V_end_pos", VEndPos},
                    {"V_score", VScore},
                    {"J_hit", JHit},
                    {"J_start_pos", JStartPos},
                    {"J_end_pos", JEndPos},
                    {"J_score", JScore},
            };
}
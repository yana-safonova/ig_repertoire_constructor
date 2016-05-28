#include "cdr_config.hpp"
#include <config_common.hpp>

namespace cdr_labeler {

    void load(CDRLabelerConfig::InputParams &ip, boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(ip.input_reads, pt, "input_reads");
        load(ip.vj_finder_config, pt, "vj_finder_config");
    }

    void load(CDRLabelerConfig::OutputParams &op, boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(op.output_dir, pt, "output_dir");
    }

    void load(CDRLabelerConfig::RunParams &rp, boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(rp.num_threads, pt, "num_threads");
    }

    void load(CDRLabelerConfig::CDRsParams::HCDR1Params &cdr1_p, boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(cdr1_p.max_length, pt, "max_length");
        load(cdr1_p.min_length, pt, "min_length");
        load(cdr1_p.residues_after, pt, "residues_after");
        load(cdr1_p.residues_before, pt, "residues_before");
        load(cdr1_p.start_pos, pt, "start_pos");
        load(cdr1_p.start_shift, pt, "start_shift");
    }

    void load(CDRLabelerConfig::CDRsParams::HCDR2Params &cdr2_p, boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(cdr2_p.max_length, pt, "max_length");
        load(cdr2_p.min_length, pt, "min_length");
        load(cdr2_p.residues_after, pt, "residues_after");
        load(cdr2_p.residues_before, pt, "residues_before");
        load(cdr2_p.distance_from_cdr1_end, pt, "distance_from_cdr1_end");
        load(cdr2_p.distance_shift, pt, "distance_shift");
    }

    void load(CDRLabelerConfig::CDRsParams::HCDR3Params &cdr3_p, boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(cdr3_p.max_length, pt, "max_length");
        load(cdr3_p.min_length, pt, "min_length");
        load(cdr3_p.residues_after, pt, "residues_after");
        load(cdr3_p.residues_before, pt, "residues_before");
        load(cdr3_p.distance_from_cdr2_end, pt, "distance_from_cdr2_end");
        load(cdr3_p.distance_shift, pt, "distance_shift");
    }

    void load(CDRLabelerConfig::CDRsParams &cdrs_p, boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(cdrs_p.hcdr1_params, pt, "hcdr1_params");
        load(cdrs_p.hcdr2_params, pt, "hcdr2_params");
        load(cdrs_p.hcdr3_params, pt, "hcdr3_params");
    }

    void CDRLabelerConfig::load(std::string config_fname) {
        boost::property_tree::ptree pt;
        boost::property_tree::read_info(config_fname, pt);
        using config_common::load;
        load(input_params, pt, "input_params", true);
        load(output_params, pt, "output_params", true);
        load(run_params, pt, "run_params", true);
        load(cdrs_params, pt, "cdrs_params", true);
        vj_finder::load(vj_finder_config, input_params.vj_finder_config);
        vj_finder_config.algorithm_params.germline_params.pseudogenes = false;
        vj_finder_config.algorithm_params.germline_params.loci = "IGH";
    }
}
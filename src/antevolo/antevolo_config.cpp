#include "antevolo_config.hpp"
#include <config_common.hpp>

namespace antevolo {
    void load(AntEvoloConfig::InputParams &input_params, boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(input_params.cdr_labeler_config_fname, pt, "cdr_labeler_config_fname");
        load(input_params.input_reads, pt, "input_reads");
        load(input_params.model_input, pt, "model_input");
    }

    void update_paths(AntEvoloConfig::OutputParams &output_params) {
        output_params.cdr_graph_dir = path::append_path(output_params.output_dir, output_params.cdr_graph_dir);
        output_params.tree_dir = path::append_path(output_params.output_dir, output_params.tree_dir);
        output_params.vertex_dir = path::append_path(output_params.output_dir, output_params.vertex_dir);
        output_params.trash_output = path::append_path(output_params.output_dir, output_params.trash_output);
    }

    void load(AntEvoloConfig::OutputParams &output_params, boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(output_params.output_dir, pt, "output_dir");
        load(output_params.cdr_graph_dir, pt, "cdr_graph_dir");
        load(output_params.trash_output, pt, "trash_output");
        load(output_params.tree_dir, pt, "tree_dir");
        load(output_params.vertex_dir, pt, "vertex_dir");
        update_paths(output_params);
    }

    void load(AntEvoloConfig::RunParams &run_params, boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(run_params.num_threads, pt, "num_threads");
    }

    void load(AntEvoloConfig::AlgorithmParams::SimilarCDR3Params &similar_cdr3s_params, boost::property_tree::ptree const &pt,
              bool) {
        using config_common::load;
        load(similar_cdr3s_params.num_mismatches, pt, "num_mismatches");
        load(similar_cdr3s_params.num_indels, pt, "num_indels");
    }

    void load(AntEvoloConfig::AlgorithmParams::EdgeConstructionParams &edge_construction_params,
              boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(edge_construction_params.min_num_intersected_v_shms, pt, "min_num_intersected_v_shms");
        load(edge_construction_params.intersected_edge_coeff, pt, "intersected_edge_coeff");
    }

    void load(AntEvoloConfig::AlgorithmParams &algorithm_params, boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(algorithm_params.similar_cdr3s_params, pt, "similar_cdr3s_params");
        load(algorithm_params.edge_construction_params, pt, "edge_construction_params");
    }

    void AntEvoloConfig::load(std::string antevolo_config_fname) {
        boost::property_tree::ptree pt;
        boost::property_tree::read_info(antevolo_config_fname, pt);
        using config_common::load;
        load(input_params, pt, "input_params");
        load(output_params, pt, "output_params");
        load(run_params, pt, "run_params");
        load(algorithm_params, pt, "algorithm_params");
        cdr_labeler_config.load(input_params.cdr_labeler_config_fname);
        cdr_labeler_config.vj_finder_config.algorithm_params.germline_params.loci = "IG";
        cdr_labeler_config.shm_params.shm_finding_algorithm =
                cdr_labeler::CDRLabelerConfig::SHMFindingParams::SHMFindingAlgorithm::CDRFilteringSHMAlgorithm;
    }
}
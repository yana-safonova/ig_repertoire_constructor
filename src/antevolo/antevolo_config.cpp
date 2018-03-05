#include "antevolo_config.hpp"
#include "shm_kmer_matrix_estimator_config.hpp"
#include <config_common.hpp>

namespace antevolo {
    void load(AntEvoloConfig::InputParams &input_params, boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(input_params.input_reads, pt, "input_reads");
        load(input_params.decomposition_rcm, pt, "decomposition_rcm");
        load(input_params.cdr_labeler_config_fname, pt, "cdr_labeler_config_fname");
        load(input_params.shm_kmer_matrix_estimator_config_fname, pt, "shm_kmer_matrix_estimator_config_fname");

        load(input_params.shm_kmer_model_igh, pt, "shm_kmer_model_igh");
        load(input_params.shm_kmer_model_igk, pt, "shm_kmer_model_igk");
        load(input_params.shm_kmer_model_igl, pt, "shm_kmer_model_igl");
    }

    void update_paths(AntEvoloConfig::OutputParams &output_params) {
        output_params.cdr_graph_dir = path::append_path(output_params.output_dir, output_params.cdr_graph_dir);
        output_params.tree_dir = path::append_path(output_params.output_dir, output_params.tree_dir);
        output_params.vertex_dir = path::append_path(output_params.output_dir, output_params.vertex_dir);
        output_params.trash_output = path::append_path(output_params.output_dir, output_params.trash_output);
        output_params.tree_details = path::append_path(output_params.output_dir, output_params.tree_details);
        output_params.tree_shm_dir = path::append_path(output_params.output_dir, output_params.tree_shm_dir);
        output_params.parallel_shm_output.parallel_bulges_dir = path::append_path(output_params.output_dir,
                                                                                  output_params.parallel_shm_output.parallel_bulges_dir);
        output_params.parallel_shm_output.parallel_shm_dir = path::append_path(output_params.output_dir,
                                                                               output_params.parallel_shm_output.parallel_shm_dir);
    }

    void load(AntEvoloConfig::OutputParams::ParallelSHMOutput &parallel_shm_output,
              boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(parallel_shm_output.parallel_bulges_dir, pt, "parallel_bulges_dir");
        load(parallel_shm_output.parallel_shm_dir, pt, "parallel_shm_dir");
    }

    void load(AntEvoloConfig::OutputParams &output_params, boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(output_params.output_dir, pt, "output_dir");
        load(output_params.cdr_graph_dir, pt, "cdr_graph_dir");
        load(output_params.trash_output, pt, "trash_output");
        load(output_params.tree_dir, pt, "tree_dir");
        load(output_params.vertex_dir, pt, "vertex_dir");
        load(output_params.tree_details, pt, "tree_details");
        load(output_params.tree_shm_dir, pt, "tree_shm_dir");
        load(output_params.parallel_shm_output, pt, "parallel_shm_output");
        update_paths(output_params);
    }

    void load(AntEvoloConfig::RunParams &run_params, boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(run_params.num_threads, pt, "num_threads");
    }

    void load(AntEvoloConfig::AlgorithmParams::SimilarCDR3Params &similar_cdr3s_params,
              boost::property_tree::ptree const &pt,
              bool) {
        using config_common::load;
        load(similar_cdr3s_params.num_mismatches_igh, pt, "num_mismatches_igh");
        load(similar_cdr3s_params.num_mismatches_igk, pt, "num_mismatches_igk");
        load(similar_cdr3s_params.num_mismatches_igl, pt, "num_mismatches_igl");
        load(similar_cdr3s_params.num_indels, pt, "num_indels");
    }

    void load(AntEvoloConfig::AlgorithmParams::EdgeConstructionParams &edge_construction_params,
              boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(edge_construction_params.min_num_intersected_v_shms, pt, "min_num_intersected_v_shms");
        load(edge_construction_params.intersected_edge_coeff, pt, "intersected_edge_coeff");
    }

    void load(AntEvoloConfig::AlgorithmParams::ParallelEvolutionParams &parallel_params,
              boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(parallel_params.enable_parallel_shms_finder, pt, "enable_parallel_shms_finder");
    }

    void load(AntEvoloConfig::AlgorithmParams &algorithm_params, boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(algorithm_params.compare, pt, "compare");
        load(algorithm_params.model, pt, "model");
        load(algorithm_params.similar_cdr3s_params, pt, "similar_cdr3s_params");
        load(algorithm_params.edge_construction_params, pt, "edge_construction_params");
        load(algorithm_params.parallel_evolution_params, pt, "parallel_evolution_params");
    }

    void AntEvoloConfig::load(std::string antevolo_config_fname) {
        boost::property_tree::ptree pt;
        boost::property_tree::read_info(antevolo_config_fname, pt);
        using config_common::load;
        load(input_params, pt, "input_params");
        load(output_params, pt, "output_params");
        load(run_params, pt, "run_params");
        load(algorithm_params, pt, "algorithm_params");
        cdr_labeler::load(cdr_labeler_config, input_params.cdr_labeler_config_fname);
        //INFO("loci:" << cdr_labeler_config.vj_finder_config.algorithm_params.germline_params.loci);
        //cdr_labeler_config.vj_finder_config.algorithm_params.germline_params.loci = "IGH";
        cdr_labeler_config.shm_params.shm_finding_algorithm =
                cdr_labeler::CDRLabelerConfig::SHMFindingParams::SHMFindingAlgorithm::CDRFilteringSHMAlgorithm;
        shm_kmer_matrix_estimator::load(shm_config, input_params.shm_kmer_matrix_estimator_config_fname);
    }

    size_t
    AntEvoloConfig::AlgorithmParams::GetNumMismatchesByChainType(const germline_utils::ImmuneChainType chain) const {
        size_t num_mismatches = 0;
        if (chain == germline_utils::ImmuneChainType::HeavyIgChain)
            num_mismatches = similar_cdr3s_params.num_mismatches_igh;
        else if (chain == germline_utils::ImmuneChainType::KappaIgChain)
            num_mismatches = similar_cdr3s_params.num_mismatches_igk;
        else if (chain == germline_utils::ImmuneChainType::LambdaIgChain)
            num_mismatches = similar_cdr3s_params.num_mismatches_igl;
        VERIFY(num_mismatches != 0);
        return num_mismatches;
    }
}
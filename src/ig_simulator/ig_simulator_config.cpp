//
// Created by Andrew Bzikadze on 3/15/17.
//

#include "ig_simulator_config.hpp"
#include <boost/property_tree/ptree_fwd.hpp>
#include <config_common.hpp>

namespace ig_simulator {

// IOParams start
void load(IgSimulatorConfig::IOParams::InputParams &input_params, boost::property_tree::ptree const &pt, bool) {
    using config_common::load;
    load(input_params.germline_input, pt, "germline_input");
    load(input_params.cdr_labeler_config_filename, pt, "cdr_labeler_config_filename");
}

void load(IgSimulatorConfig::IOParams::OutputParams &output_params, boost::property_tree::ptree const &pt, bool) {
    using config_common::load;
    load(output_params.log_filename, pt, "log_filename");
    load(output_params.output_dir, pt, "output_dir");
}

void load(IgSimulatorConfig::IOParams &io_params, boost::property_tree::ptree const &pt, bool) {
    using config_common::load;
    load(io_params.input_params, pt, "input_params");
    load(io_params.output_params, pt, "output_params");
}
// IOParams end

// AlgorithmParams start
void load(IgSimulatorConfig::AlgorithmParams &algorithm_params,
          boost::property_tree::ptree const &pt, bool)
{
    using config_common::load;
    load(algorithm_params.germline_params, pt, "germline_params");
}
// AlgorithmParams end

// SimulationParams start
void load(IgSimulatorConfig::SimulationParams::BaseRepertoireParams::MetarootSimulationParams
          ::GeneChooserParams &gene_chooser_params,
          boost::property_tree::ptree const &pt, bool) {
    using config_common::load;
    using GeneChooserMethod =
        IgSimulatorConfig::SimulationParams::BaseRepertoireParams::
        MetarootSimulationParams::GeneChooserParams::GeneChooserMethod;
    std::string method_str(pt.get<std::string>("gene_chooser_method"));
    std::string method_str_lowercase(method_str);
    std::transform(method_str.begin(), method_str.end(),
                   method_str_lowercase.begin(), ::tolower);
    if (method_str == "uniform") {
        gene_chooser_params.method = GeneChooserMethod::Uniform;
    }
}

void load(IgSimulatorConfig::SimulationParams::BaseRepertoireParams::MetarootSimulationParams
          ::NucleotidesRemoverParams::UniformRemoverParams &uniform_remover_params,
          boost::property_tree::ptree const &pt, bool) {
    using config_common::load;
    load(uniform_remover_params.max_remove_v_gene, pt, "max_remove_v_gene");
    load(uniform_remover_params.max_remove_d_gene_left, pt, "max_remove_d_gene_left");
    load(uniform_remover_params.max_remove_d_gene_right, pt, "max_remove_d_gene_right");
    load(uniform_remover_params.max_remove_j_gene, pt, "max_remove_j_gene");
}

void load(IgSimulatorConfig::SimulationParams::BaseRepertoireParams::MetarootSimulationParams
          ::NucleotidesRemoverParams &nucleotides_remover_params,
          boost::property_tree::ptree const &pt, bool) {
    using config_common::load;
    using NucleotidesRemoverMethod =
        IgSimulatorConfig::SimulationParams::BaseRepertoireParams::
        MetarootSimulationParams::NucleotidesRemoverParams::NucleotidesRemoverMethod;
    std::string method_str(pt.get<std::string>("nucleotides_remover_method"));
    std::string method_str_lowercase(method_str);
    std::transform(method_str.begin(), method_str.end(),
                   method_str_lowercase.begin(), ::tolower);
    if (method_str == "uniform") {
        nucleotides_remover_params.method = NucleotidesRemoverMethod::Uniform;
        load(nucleotides_remover_params.uniform_remover_params, pt, "uniform_remover_params");
    }
}

void load(IgSimulatorConfig::SimulationParams::BaseRepertoireParams::MetarootSimulationParams
          ::PNucleotidesCreatorParams::UniformCreatorParams &uniform_creator_params,
          boost::property_tree::ptree const &pt, bool) {
    using config_common::load;
    load(uniform_creator_params.max_create_v_gene, pt, "max_create_v_gene");
    load(uniform_creator_params.max_create_d_gene_left, pt, "max_create_d_gene_left");
    load(uniform_creator_params.max_create_d_gene_right, pt, "max_create_d_gene_right");
    load(uniform_creator_params.max_create_j_gene, pt, "max_create_j_gene");
}

void load(IgSimulatorConfig::SimulationParams::BaseRepertoireParams::MetarootSimulationParams
          ::PNucleotidesCreatorParams &p_nucleptides_creator_params,
          boost::property_tree::ptree const &pt, bool) {
    using config_common::load;
    using PNucleotidesCreatorParams =
        IgSimulatorConfig::SimulationParams::BaseRepertoireParams::
        MetarootSimulationParams::PNucleotidesCreatorParams::PNucleotidesCreatorMethod;
    std::string method_str(pt.get<std::string>("p_nucleotides_creator_method"));
    std::string method_str_lowercase(method_str);
    std::transform(method_str.begin(), method_str.end(),
                   method_str_lowercase.begin(), ::tolower);
    if (method_str == "uniform") {
        p_nucleptides_creator_params.method = PNucleotidesCreatorParams::Uniform;
        load(p_nucleptides_creator_params.uniform_creator_params, pt, "uniform_creator_params");
    }
}

void load(IgSimulatorConfig::SimulationParams::BaseRepertoireParams::MetarootSimulationParams
          ::NNucleotidesInserterParams::UniformInserterParams &uniform_inserter_params,
          boost::property_tree::ptree const &pt, bool) {
    using config_common::load;
    load(uniform_inserter_params.max_vj_insertion, pt, "max_vj_insertion");
    load(uniform_inserter_params.max_vd_insertion, pt, "max_vd_insertion");
    load(uniform_inserter_params.max_dj_insertion, pt, "max_dj_insertion");
}

void load(IgSimulatorConfig::SimulationParams::BaseRepertoireParams::MetarootSimulationParams
          ::NNucleotidesInserterParams &n_nucleotides_inserter_params,
          boost::property_tree::ptree const &pt, bool) {
    using config_common::load;
    using NNucleotidesInserterParams =
        IgSimulatorConfig::SimulationParams::BaseRepertoireParams::
        MetarootSimulationParams::NNucleotidesInserterParams::NNucleotidesInserterMethod;
    std::string method_str(pt.get<std::string>("n_nucleotides_method"));
    std::string method_str_lowercase(method_str);
    std::transform(method_str.begin(), method_str.end(),
                   method_str_lowercase.begin(), ::tolower);
    if (method_str == "uniform") {
        n_nucleotides_inserter_params.method = NNucleotidesInserterParams::Uniform;
        load(n_nucleotides_inserter_params.uniform_inserter_params, pt, "uniform_inserter_params");
    }
}

void load(IgSimulatorConfig::SimulationParams::BaseRepertoireParams::MetarootSimulationParams
          ::CleavageParams &cleavage_params,
          boost::property_tree::ptree const &pt, bool) {
    using config_common::load;
    load(cleavage_params.prob_cleavage_v, pt, "prob_cleavage_v");
    load(cleavage_params.prob_cleavage_d_left, pt, "prob_cleavage_d_left");
    load(cleavage_params.prob_cleavage_d_right, pt, "prob_cleavage_d_right");
    load(cleavage_params.prob_cleavage_j, pt, "prob_cleavage_j");
}

void load(IgSimulatorConfig::SimulationParams::BaseRepertoireParams
          ::MetarootSimulationParams &metaroot_simulation_params,
          boost::property_tree::ptree const &pt, bool)
{
    using config_common::load;
    load(metaroot_simulation_params.gene_chooser_params, pt, "gene_chooser_params");
    load(metaroot_simulation_params.nucleotides_remover_params, pt, "nucleotides_remover_params");
    load(metaroot_simulation_params.p_nucleotides_creator_params, pt, "p_nucleotides_creator_params");
    load(metaroot_simulation_params.n_nucleotides_inserter_params, pt, "n_nucleotides_inserter_params");
    load(metaroot_simulation_params.cleavage_params, pt, "cleavage_params");
}

void load(IgSimulatorConfig::SimulationParams::BaseRepertoireParams &base_repertoire_params,
          boost::property_tree::ptree const &pt, bool)
{
    using config_common::load;
    load(base_repertoire_params.metaroot_simulation_params, pt, "metaroot_simulation_params");
}

void load(IgSimulatorConfig::SimulationParams &simulation_params,
          boost::property_tree::ptree const &pt, bool)
{
    using config_common::load;
    load(simulation_params.base_repertoire_params, pt, "base_repertoire_params");
}
// SimulationParams end


void load(IgSimulatorConfig &cfg, boost::property_tree::ptree const &pt, bool complete) {
    using config_common::load;
    load(cfg.io_params, pt, "io_params", complete);
    load(cfg.algorithm_params, pt, "algorithm_params", complete);
    load(cfg.simulation_params, pt, "simulation_params", complete);
    // TODO remove this hack
    cfg.simulation_params.base_repertoire_params.
        metaroot_simulation_params.cdr_labeler_config.load(cfg.io_params.input_params.cdr_labeler_config_filename);
}

void load(IgSimulatorConfig &cfg, std::string const &filename) {
    boost::property_tree::ptree pt;
    boost::property_tree::read_info(filename, pt);
    load(cfg, pt, true);
}

} // End namespace ig_simulator

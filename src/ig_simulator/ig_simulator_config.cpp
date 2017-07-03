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
    load(output_params.base_repertoire_filename, pt, "base_repertoire_filename");
    load(output_params.base_repertoire_info, pt, "base_repertoire_info");
    load(output_params.filtered_pool, pt, "filtered_pool");
    load(output_params.full_pool, pt, "full_pool");
    load(output_params.trees_dir, pt, "trees_dir");
}

void load(IgSimulatorConfig::IOParams &io_params, boost::property_tree::ptree const &pt, bool) {
    using config_common::load;
    load(io_params.input_params, pt, "input_params");
    load(io_params.output_params, pt, "output_params");
}
// IOParams end

// SimulationParams start
void load(GeneChooserParams::CustomGeneChooserParams& custom_gene_chooser_params,
          boost::property_tree::ptree const &pt, bool)
{
    using config_common::load;
    load(custom_gene_chooser_params.v_genes_probs, pt, "v_genes_probs");
    load(custom_gene_chooser_params.v_genes_probs, pt, "d_genes_probs");
    load(custom_gene_chooser_params.v_genes_probs, pt, "j_genes_probs");
}


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
    } else if (method_str == "custom") {
        gene_chooser_params.method = GeneChooserMethod::Custom;
        load(gene_chooser_params.custom_gene_chooser_params, pt, "custom_chooser_params");
    } else {
        VERIFY(false);
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
    } else {
        VERIFY(false);
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
    } else {
        VERIFY(false);
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
    } else {
        VERIFY(false);
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

void load(MultiplicityCreatorParams::GeometricParams &geometric_params,
          boost::property_tree::ptree const &pt, bool)
{
    using config_common::load;
    load(geometric_params.lambda, pt, "lambda");
}

void load(MultiplicityCreatorParams &multiplicity_creator_params,
          boost::property_tree::ptree const &pt, bool)
{
    using config_common::load;
    using MultiplicityCreatorMethod = MultiplicityCreatorParams::MultiplicityCreatorMethod;

    std::string method_str(pt.get<std::string>("multiplicity_method"));
    std::string method_str_lowercase(method_str);
    std::transform(method_str.begin(), method_str.end(),
                   method_str_lowercase.begin(), ::tolower);
    if (method_str == "geometric") {
        multiplicity_creator_params.method = MultiplicityCreatorMethod::Geometric;
        load(multiplicity_creator_params.geometric_params, pt, "geometric_params");
    } else {
        VERIFY(false);
    }
}

void load(ProductiveParams &base_repertoire_params,
          boost::property_tree::ptree const &pt, bool)
{
    using config_common::load;
    load(base_repertoire_params.productive_part, pt, "productive_part");
}

void load(IgSimulatorConfig::SimulationParams::BaseRepertoireParams &base_repertoire_params,
          boost::property_tree::ptree const &pt, bool)
{
    using config_common::load;
    load(base_repertoire_params.metaroot_simulation_params, pt, "metaroot_simulation_params");
    load(base_repertoire_params.multiplicity_creator_params, pt, "multiplicity_creator_params");
    load(base_repertoire_params.productive_params, pt, "productive_params");
    load(base_repertoire_params.number_of_metaroots, pt, "number_of_metaroots");
}

void load(TreeSizeGeneratorParams::GeometricParams &geometric_params,
          boost::property_tree::ptree const &pt, bool)
{
    using config_common::load;
    load(geometric_params.lambda, pt, "lambda");
}

void load(TreeSizeGeneratorParams &tree_size_generator_params,
          boost::property_tree::ptree const &pt, bool)
{
    using config_common::load;
    using TreeSizeGeneratorMethod = TreeSizeGeneratorParams::TreeSizeGeneratorMethod;

    std::string method_str(pt.get<std::string>("tree_size_generator_method"));
    std::string method_str_lowercase(method_str);
    std::transform(method_str.begin(), method_str.end(),
                   method_str_lowercase.begin(), ::tolower);
    if (method_str == "geometric") {
        tree_size_generator_params.method = TreeSizeGeneratorMethod::Geometric;
        load(tree_size_generator_params.geometric_params, pt, "geometric_params");
    } else {
        VERIFY(false);
    }
}

void load(SHM_CreatorParams::PoissonCreatorParams &params,
          boost::property_tree::ptree const &pt, bool)
{
    using config_common::load;
    load(params.lambda, pt, "lambda");
}

void load(SHM_CreatorParams &shm_creator_params,
          boost::property_tree::ptree const &pt, bool) {
    using config_common::load;
    using SHM_CreatorMethod = SHM_CreatorParams::SHM_CreatorMethod;

    std::string method_str(pt.get<std::string>("shm_creator_method"));
    std::string method_str_lowercase(method_str);
    std::transform(method_str.begin(), method_str.end(),
                   method_str_lowercase.begin(), ::tolower);
    if (method_str == "poisson") {
        shm_creator_params.method = SHM_CreatorMethod::Poisson;
        load(shm_creator_params.poisson_params, pt, "poisson_params");
    } else {
        VERIFY(false);
    }
}

void load(ClonalTreeSimulatorParams &clonal_tree_simulator_params,
          boost::property_tree::ptree const &pt, bool)
{
    using config_common::load;
    using PoolManagerStrategy = ClonalTreeSimulatorParams::PoolManagerStrategy;

    std::string method_str(pt.get<std::string>("pool_manager_strategy"));
    std::string method_str_lowercase(method_str);
    std::transform(method_str.begin(), method_str.end(),
                   method_str_lowercase.begin(), ::tolower);
    if (method_str == "uniform") {
        clonal_tree_simulator_params.pool_manager_strategy = PoolManagerStrategy::UniformPoolManager;
    } else if (method_str == "wide") {
        clonal_tree_simulator_params.pool_manager_strategy = PoolManagerStrategy::WideTreePoolManager;
    } else if (method_str == "deep") {
        clonal_tree_simulator_params.pool_manager_strategy = PoolManagerStrategy::DeepTreePoolManager;
    } else {
        VERIFY(false);
    }

    load(clonal_tree_simulator_params.prob_ret_to_pool, pt, "prob_ret_to_pool");
    load(clonal_tree_simulator_params.lambda_distr_n_children, pt, "lambda_distr_n_children");
    load(clonal_tree_simulator_params.tree_size_generator_params, pt, "tree_size_generator_params");
    load(clonal_tree_simulator_params.shm_creator_params, pt, "shm_creator_params");
}

void load(IgSimulatorConfig::SimulationParams &simulation_params,
          boost::property_tree::ptree const &pt, bool)
{
    using config_common::load;
    load(simulation_params.base_repertoire_params, pt, "base_repertoire_params");
    load(simulation_params.clonal_tree_simulator_params, pt, "clonal_tree_simulator_params");
}
// SimulationParams end


void load(IgSimulatorConfig &cfg, boost::property_tree::ptree const &pt, bool complete) {
    using config_common::load;
    load(cfg.io_params, pt, "io_params", complete);
    load(cfg.simulation_params, pt, "simulation_params", complete);
    load(cfg.germline_params, pt, "germline_params");
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

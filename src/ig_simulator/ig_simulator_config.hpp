//
// Created by Andrew Bzikadze on 3/15/17.
//

#pragma once

#include "io/library.hpp"
#include <logger/logger.hpp>
#include "config_singl.hpp"
#include "germline_utils/germline_config.hpp"
#include "cdr_config.hpp"

namespace ig_simulator {

struct IgSimulatorConfig {
    struct IOParams {
        struct InputParams {
            germline_utils::GermlineInput germline_input;
            std::string cdr_labeler_config_filename;
        };

        struct OutputParams {
            std::string output_dir;
            std::string log_filename;
            std::string base_repertoire_filename;
            std::string base_repertoire_info;
            std::string filtered_pool;
            std::string full_pool;
            std::string trees_dir;
        };

        InputParams input_params;
        OutputParams output_params;
    };


    struct SimulationParams {
        struct BaseRepertoireParams {
            struct MetarootSimulationParams {
                struct GeneChooserParams {
                    struct CustomGeneChooserParams {
                        std::string v_genes_probs;
                        std::string d_genes_probs;
                        std::string j_genes_probs;
                    };

                    enum class GeneChooserMethod { Uniform, Custom };
                    GeneChooserMethod method;
                    CustomGeneChooserParams custom_gene_chooser_params;
                };

                struct NucleotidesRemoverParams {
                    enum class NucleotidesRemoverMethod { Uniform };
                    struct UniformRemoverParams {
                        size_t max_remove_v_gene;
                        size_t max_remove_d_gene_left;
                        size_t max_remove_d_gene_right;
                        size_t max_remove_j_gene;
                    };
                    NucleotidesRemoverMethod method;
                    UniformRemoverParams uniform_remover_params;
                };

                struct PNucleotidesCreatorParams {
                    enum class PNucleotidesCreatorMethod { Uniform };
                    struct UniformCreatorParams {
                        size_t max_create_v_gene;
                        size_t max_create_d_gene_left;
                        size_t max_create_d_gene_right;
                        size_t max_create_j_gene;
                    };
                    PNucleotidesCreatorMethod method;
                    UniformCreatorParams uniform_creator_params;
                };

                struct NNucleotidesInserterParams {
                    enum class NNucleotidesInserterMethod { Uniform };
                    struct UniformInserterParams {
                        size_t max_vj_insertion;
                        size_t max_vd_insertion;
                        size_t max_dj_insertion;
                    };
                    NNucleotidesInserterMethod method;
                    UniformInserterParams uniform_inserter_params;
                };

                struct CleavageParams {
                    double prob_cleavage_v;
                    double prob_cleavage_d_left;
                    double prob_cleavage_d_right;
                    double prob_cleavage_j;
                };

                GeneChooserParams gene_chooser_params;
                NucleotidesRemoverParams nucleotides_remover_params;
                PNucleotidesCreatorParams p_nucleotides_creator_params;
                NNucleotidesInserterParams n_nucleotides_inserter_params;
                CleavageParams cleavage_params;
                cdr_labeler::CDRLabelerConfig cdr_labeler_config;
            };

            struct MultiplicityCreatorParams {
                struct GeometricParams {
                    double lambda;
                };

                enum class MultiplicityCreatorMethod { Geometric };
                MultiplicityCreatorMethod method;
                GeometricParams geometric_params;
            };

            struct ProductiveParams {
                double productive_part;
            };

            MetarootSimulationParams metaroot_simulation_params;
            MultiplicityCreatorParams multiplicity_creator_params;
            ProductiveParams productive_params;

            size_t number_of_metaroots;
        };

        struct ClonalTreeSimulatorParams {
            struct TreeSizeGeneratorParams {
                struct GeometricParams {
                    double lambda;
                };

                enum class TreeSizeGeneratorMethod { Geometric };
                TreeSizeGeneratorMethod method;
                GeometricParams geometric_params;
            };

            struct SHM_CreatorParams {
                struct PoissonCreatorParams {
                    double lambda;
                };

                enum class SHM_CreatorMethod { Poisson };
                SHM_CreatorMethod method;
                PoissonCreatorParams poisson_params;
            };

            enum class PoolManagerStrategy { UniformPoolManager, WideTreePoolManager, DeepTreePoolManager };
            PoolManagerStrategy pool_manager_strategy;

            double prob_ret_to_pool;
            double lambda_distr_n_children;
            TreeSizeGeneratorParams tree_size_generator_params;
            SHM_CreatorParams shm_creator_params;
        };

        BaseRepertoireParams base_repertoire_params;
        ClonalTreeSimulatorParams clonal_tree_simulator_params;
    };

    IOParams io_params;
    germline_utils::GermlineParams germline_params;
    SimulationParams simulation_params;
};

using BaseRepertoireParams = IgSimulatorConfig::SimulationParams::BaseRepertoireParams;
using ClonalTreeSimulatorParams = IgSimulatorConfig::SimulationParams::ClonalTreeSimulatorParams;

using MetarootSimulationParams = BaseRepertoireParams::MetarootSimulationParams;
using MultiplicityCreatorParams = BaseRepertoireParams::MultiplicityCreatorParams;
using ProductiveParams = BaseRepertoireParams::ProductiveParams;

using MultiplicityCreatorMethod = MultiplicityCreatorParams::MultiplicityCreatorMethod;

using GeneChooserParams = MetarootSimulationParams::GeneChooserParams;
using GeneChooserMethod = GeneChooserParams::GeneChooserMethod;

using NucleotidesRemoverParams = MetarootSimulationParams::NucleotidesRemoverParams;
using NucleotidesRemoverMethod = NucleotidesRemoverParams::NucleotidesRemoverMethod;

using PNucleotidesCreatorParams = MetarootSimulationParams::PNucleotidesCreatorParams;
using PNucleotidesCreatorMethod = PNucleotidesCreatorParams::PNucleotidesCreatorMethod;

using NNucleotidesInserterParams = MetarootSimulationParams::NNucleotidesInserterParams;
using NNucleotidesInserterMethod = NNucleotidesInserterParams::NNucleotidesInserterMethod;

using CleavageParams = MetarootSimulationParams::CleavageParams;

using TreeSizeGeneratorParams = ClonalTreeSimulatorParams::TreeSizeGeneratorParams;
using SHM_CreatorParams = ClonalTreeSimulatorParams::SHM_CreatorParams;
using PoolManagerStrategy = ClonalTreeSimulatorParams::PoolManagerStrategy;

void load(IgSimulatorConfig &cfg, std::string const &filename);

typedef config_common::config<IgSimulatorConfig> igs_cfg;

} // End namespace ig_simulator

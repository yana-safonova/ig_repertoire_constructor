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
        };

        InputParams input_params;
        OutputParams output_params;
    };

    struct AlgorithmParams {
        germline_utils::GermlineParams germline_params;
    };

    struct SimulationParams {
        struct BaseRepertoireParams {
            struct MetarootSimulationParams {
                struct GeneChooserParams {
                    enum class GeneChooserMethod { Uniform };
                    GeneChooserMethod method;
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
            MetarootSimulationParams metaroot_simulation_params;
            MultiplicityCreatorParams multiplicity_creator_params;
        };

        BaseRepertoireParams base_repertoire_params;
    };

    IOParams io_params;
    AlgorithmParams algorithm_params;
    SimulationParams simulation_params;
};

using MetarootSimulationParams = IgSimulatorConfig::SimulationParams::BaseRepertoireParams::MetarootSimulationParams;
using MultiplicityCreatorParams = IgSimulatorConfig::SimulationParams::BaseRepertoireParams::MultiplicityCreatorParams;

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


void load(IgSimulatorConfig &cfg, std::string const &filename);

typedef config_common::config<IgSimulatorConfig> igs_cfg;

} // End namespace ig_simulator

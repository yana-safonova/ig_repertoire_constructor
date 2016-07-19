#pragma once

#include "io/library.hpp"
#include "config_singl.hpp"
#include <boost/property_tree/ptree_fwd.hpp>

#include "../vj_finder/vj_finder_config.hpp"

namespace vdj_labeler {

struct VDJLabelerConfig {
    struct IOParams {
        struct InputParams {
            std::string input_sequences;
            std::string vj_finder_config;

            struct GermlineGenes {
                struct IGH_Genes {
                    std::string variable_genes;
                    std::string diversity_genes;
                    std::string join_genes;
                };

                IGH_Genes igh_genes;
            };

            GermlineGenes germlines;
        };

        struct OutputParams {
            std::string log_filename;
            std::string output_dir;
        };

        InputParams input_params;
        OutputParams output_params;
    };

    struct RunParams {
        unsigned threads_count;
        unsigned max_memory;
    };

    struct DAlignmentQualityParams {
        unsigned min_coverage;
    };

    IOParams io_params;
    RunParams run_params;
    DAlignmentQualityParams d_align_quality_params;
    vj_finder::VJFinderConfig vj_finder_config;

    void load(const std::string &filename);
};

} // End namespace vdj_labeler

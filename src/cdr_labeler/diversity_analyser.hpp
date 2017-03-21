#pragma once

#include <../graph_utils/sparse_graph.hpp>
#include <annotation_utils/annotated_clone_set.hpp>
#include "compressed_cdr_set.hpp"
#include "cdr_config.hpp"

namespace cdr_labeler {
    class DiversityAnalyser {
        //const annotation_utils::CDRAnnotatedCloneSet &clone_set_;
        const CDRLabelerConfig::InputParams &input_params_;
        const CDRLabelerConfig::OutputParams &output_params_;

        CompressedCDRSet cdr1_compressed_set_;
        CompressedCDRSet cdr2_compressed_set_;
        CompressedCDRSet cdr3_compressed_set_;
        std::vector<SparseGraphPtr> cdr3_graphs_;
        GraphComponentMap graph_component_map_;

        void InitializeGraph(std::string compressed_cdr3_fasta);

        //size_t ComputeD50(const std::vector<SparseGraphPtr> connected_components) const;

        const CompressedCDRSet& GetCompressedCloneSet(annotation_utils::StructuralRegion);

        size_t MaxConnectedComponentAbundance();

    public:
        DiversityAnalyser(const annotation_utils::CDRAnnotatedCloneSet &clone_set,
                          const CDRLabelerConfig::InputParams &input_params,
                          const CDRLabelerConfig::OutputParams &output_params,
                          std::string compressed_cdr3_fasta = "") :
                //clone_set_(clone_set),
                input_params_(input_params),
                output_params_(output_params),
                cdr1_compressed_set_(annotation_utils::StructuralRegion::CDR1, clone_set),
                cdr2_compressed_set_(annotation_utils::StructuralRegion::CDR2, clone_set),
                cdr3_compressed_set_(annotation_utils::StructuralRegion::CDR3, clone_set){
            InitializeGraph(compressed_cdr3_fasta);
        }

        double ShannonIndex(annotation_utils::StructuralRegion region);

        double SimpsonIndex(annotation_utils::StructuralRegion region);

        double ClonalShannonIndex();

        double ClonalSimpsonIndex();
    };
}
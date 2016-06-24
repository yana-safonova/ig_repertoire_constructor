#pragma once

#include <../graph_utils/sparse_graph.hpp>
#include <annotation_utils/annotated_clone_set.hpp>
#include "compressed_cdr_set.hpp"

namespace cdr_labeler {
    class DiversityAnalyser {
        const annotation_utils::CDRAnnotatedCloneSet &clone_set_;

        CompressedCDRSet cdr1_compressed_set_;
        CompressedCDRSet cdr2_compressed_set_;
        CompressedCDRSet cdr3_compressed_set_;
        //SparseGraphPtr cdr_graph_;

        //void InitializeGraph(std::string compressed_cdr3_fasta);

        //size_t ComputeD50(const std::vector<SparseGraphPtr> connected_components) const;

        const CompressedCDRSet& GetCompressedCloneSet(annotation_utils::StructuralRegion);

    public:
        DiversityAnalyser(const annotation_utils::CDRAnnotatedCloneSet &clone_set) :
                clone_set_(clone_set),
                cdr1_compressed_set_(annotation_utils::StructuralRegion::CDR1, clone_set),
                cdr2_compressed_set_(annotation_utils::StructuralRegion::CDR2, clone_set),
                cdr3_compressed_set_(annotation_utils::StructuralRegion::CDR3, clone_set){
        }

        double ShannonIndex(annotation_utils::StructuralRegion region);

        double SimpsonIndex(annotation_utils::StructuralRegion region);
    };
}
#pragma once

#include "cdr_config.hpp"
#include <annotation_utils/annotated_clone_set.hpp>
#include <vj_alignment_info.hpp>
#include "compressed_cdr_set.hpp"

namespace cdr_labeler {
    class CDRLabelingWriter {
        const CDRLabelerConfig::OutputParams &output_config_;
        //const vj_finder::VJAlignmentInfo &alignment_info_;
        const annotation_utils::CDRAnnotatedCloneSet &clone_set_;

        std::ostream& OutputCloneRegion(std::ostream& out, const annotation_utils::AnnotatedClone &clone,
                                       annotation_utils::StructuralRegion region) const;

        void OutputRegionFasta(std::string output_fname, annotation_utils::StructuralRegion region) const;

        seqan::CharString GetCompressedRegionFname(annotation_utils::StructuralRegion region,
                                                   CDRKey cdr_key, size_t abundance) const;

        void OutputSHMsForRead(std::ostream& out, const annotation_utils::GeneSegmentSHMs &shms) const;

    public:
        CDRLabelingWriter(const CDRLabelerConfig::OutputParams &output_config,
                          //const vj_finder::VJAlignmentInfo &alignment_info,
                          const annotation_utils::CDRAnnotatedCloneSet &clone_set) : output_config_(output_config),
                                                                                     //alignment_info_(alignment_info),
                                                                                     clone_set_(clone_set) { }

        void OutputCDRDetails() const;

        void OutputCDR1Fasta() const;

        void OutputCDR2Fasta() const;

        void OutputCDR3Fasta() const;

        void OutputCompressedCDR3Fasta() const;

        void OutputVGeneAlignment() const;

        void OutputSHMs() const;
    };
}
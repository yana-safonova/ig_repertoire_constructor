#pragma once

#include "immune_gene_labeler.hpp"

namespace cdr_labeler {
    class AnnotatedVGeneCDRLabeler : public BaseImmuneGeneCDRLabeler {
        const CDRLabelerConfig::CDRsParams::AnnotatedSearchParams &search_params_;
        germline_utils::SegmentType segment_type_;

        struct VGeneAnnotation {
            std::string name;
            size_t cdr1_start;
            size_t cdr1_end;
            size_t cdr2_start;
            size_t cdr2_end;
            size_t cdr3_start;

            VGeneAnnotation() : name(), cdr1_start(), cdr1_end(),
                                cdr2_start(), cdr2_end(), cdr3_start() { }

            VGeneAnnotation(std::string name, size_t cdr1_start, size_t cdr1_end,
                            size_t cdr2_start, size_t cdr2_end, size_t cdr3_start) :
                    name(name), cdr1_start(cdr1_start), cdr1_end(cdr1_end),
                    cdr2_start(cdr2_start), cdr2_end(cdr2_end), cdr3_start(cdr3_start) { }
        };

        std::vector<AnnotatedVGeneCDRLabeler::VGeneAnnotation> v_annotations_;
        std::unordered_map<std::string, size_t> v_name_index_map_;

        void Initialize();

    public:
        AnnotatedVGeneCDRLabeler(const CDRLabelerConfig::CDRsParams::AnnotatedSearchParams &search_params) :
                search_params_(search_params),
                segment_type_(germline_utils::SegmentType::VariableSegment) {
            Initialize();
        }

        annotation_utils::CDRLabeling ComputeLabeling(const germline_utils::ImmuneGene &immune_gene);
    };

    class AnnotatedJGeneCDRLabeler : public BaseImmuneGeneCDRLabeler {
        const CDRLabelerConfig::CDRsParams::AnnotatedSearchParams &search_params_;
        germline_utils::SegmentType segment_type_;

        // todo: make ORF enum
        struct JGeneAnnotation {
            std::string name;
            size_t cdr3_end;
            unsigned orf;

            JGeneAnnotation() : name(), cdr3_end(), orf() { }

            JGeneAnnotation(std::string name, size_t cdr3_end) : name(name), cdr3_end(cdr3_end) {
                orf = (cdr3_end + 1) % 3;
            }
        };

        std::vector<AnnotatedJGeneCDRLabeler::JGeneAnnotation> j_annotations_;
        std::unordered_map<std::string, size_t> j_name_index_map_;

        void Initialize();

    public:
        AnnotatedJGeneCDRLabeler(const CDRLabelerConfig::CDRsParams::AnnotatedSearchParams &search_params) :
                search_params_(search_params),
                segment_type_(germline_utils::SegmentType::JoinSegment) {
            Initialize();
        }

        annotation_utils::CDRLabeling ComputeLabeling(const germline_utils::ImmuneGene &immune_gene);
    };
}
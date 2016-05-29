#pragma once

#include "hcdr1_labeler.hpp"
#include "hcdr2_labeler.hpp"
#include "hcdr3_v_labeler.hpp"
#include "hcdr3_j_labeler.hpp"

namespace cdr_labeler {
    class ImmuneGeneCDRHelper {
    public:
        static SingleLoopLabelerPtr GetCDR1Labeler(germline_utils::ImmuneGeneType gene_type,
                                                   const CDRLabelerConfig::CDRsParams &cdr_params);

        static SingleLoopLabelerPtr GetCDR2Labeler(germline_utils::ImmuneGeneType gene_type,
                                                   const CDRLabelerConfig::CDRsParams &cdr_params);

        static SingleLoopLabelerPtr GetCDR3Labeler(germline_utils::ImmuneGeneType gene_type,
                                                   const CDRLabelerConfig::CDRsParams &cdr_params);
    };

    class BaseImmuneGeneCDRLabeler {
    public:
        virtual CDRLabeling ComputeLabeling(const germline_utils::ImmuneGene &immune_gene) = 0;

        ~BaseImmuneGeneCDRLabeler() { }
    };

    class AnnotatedVGeneCDRLabeler : public BaseImmuneGeneCDRLabeler {
        const CDRLabelerConfig::CDRsParams::AnnotatedSearchParams &search_params_;
        germline_utils::SegmentType segment_type_;

        void Initialize();

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


    public:
        AnnotatedVGeneCDRLabeler(const CDRLabelerConfig::CDRsParams::AnnotatedSearchParams &search_params) :
                search_params_(search_params),
                segment_type_(germline_utils::SegmentType::VariableSegment) {
            Initialize();
        }

        CDRLabeling ComputeLabeling(const germline_utils::ImmuneGene &immune_gene);
    };

    class AnnotatedJGeneCDRLabeler : public BaseImmuneGeneCDRLabeler {
        const CDRLabelerConfig::CDRsParams::AnnotatedSearchParams &search_params_;
        germline_utils::SegmentType segment_type_;

        void Initialize();

    public:
        AnnotatedJGeneCDRLabeler(const CDRLabelerConfig::CDRsParams::AnnotatedSearchParams &search_params) :
                search_params_(search_params),
                segment_type_(germline_utils::SegmentType::JoinSegment) {
            Initialize();
        }

        CDRLabeling ComputeLabeling(const germline_utils::ImmuneGene &immune_gene);
    };

    class DeNovoImmuneGeneCDRLabeler : public BaseImmuneGeneCDRLabeler {
        germline_utils::ImmuneGeneType gene_type_;
        SingleLoopLabelerPtr cdr1_labeler_ptr_;
        SingleLoopLabelerPtr cdr2_labeler_ptr_;
        SingleLoopLabelerPtr cdr3_labeler_ptr_;

    public:
        DeNovoImmuneGeneCDRLabeler(germline_utils::ImmuneGeneType gene_type,
                                   const CDRLabelerConfig::CDRsParams &cdrs_params) :
                gene_type_(gene_type),
                cdr1_labeler_ptr_(ImmuneGeneCDRHelper::GetCDR1Labeler(gene_type, cdrs_params)),
                cdr2_labeler_ptr_(ImmuneGeneCDRHelper::GetCDR2Labeler(gene_type, cdrs_params)),
                cdr3_labeler_ptr_(ImmuneGeneCDRHelper::GetCDR3Labeler(gene_type, cdrs_params)) { }

        CDRLabeling ComputeLabeling(const germline_utils::ImmuneGene &immune_gene);
    };
}
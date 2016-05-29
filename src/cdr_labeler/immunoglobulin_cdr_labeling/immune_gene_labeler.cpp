#include "immune_gene_labeler.hpp"

#include <convert.hpp>

#include <boost/algorithm/string.hpp>

namespace cdr_labeler {
    void AnnotatedVGeneCDRLabeler::Initialize() {
        std::string v_gene_annotation = (search_params_.domain_system == CDRLabelerConfig::CDRsParams::
                                                                         AnnotatedSearchParams::DomainSystem::IMGT) ?
                            search_params_.imgt_v_annotation : search_params_.kabat_v_annotation;
        VERIFY_MSG(path::check_existence(v_gene_annotation), "File with V gene annotation " << v_gene_annotation <<
                " does not exist");
        std::ifstream ifhandler(v_gene_annotation);
        while(!ifhandler.eof()) {
            std::string tmp;
            std::getline(ifhandler, tmp);
            std::vector<std::string> splits;
            boost::split(splits, tmp, boost::is_any_of("\t"));
            AnnotatedVGeneCDRLabeler::VGeneAnnotation v_annotation(splits[search_params_.v_gene_line_index],
                            core::convert<std::string, size_t>(splits[search_params_.cdr1_start_line_index]) - 1,
                            core::convert<std::string, size_t>(splits[search_params_.cdr1_end_line_index]) - 1,
                            core::convert<std::string, size_t>(splits[search_params_.cdr2_start_line_index]) - 1,
                            core::convert<std::string, size_t>(splits[search_params_.cdr2_end_line_index]) - 1,
                            core::convert<std::string, size_t>(splits[search_params_.fr3_end_index]));
            v_annotations_.push_back(v_annotation);
            v_name_index_map_[v_annotation.name] = v_annotations_.size() - 1;
        }
    }

    CDRLabeling AnnotatedVGeneCDRLabeler::ComputeLabeling(const germline_utils::ImmuneGene &immune_gene) {
        VERIFY_MSG(immune_gene.GeneType().Segment(), "Type of immune gene " << immune_gene.name() << " is not V");
        std::string gene_name = std::string(seqan::toCString(immune_gene.name()));
        if(v_name_index_map_.find(gene_name) != v_name_index_map_.end()) {
            AnnotatedVGeneCDRLabeler::VGeneAnnotation annotation = v_annotations_[v_name_index_map_.at(gene_name)];
            return CDRLabeling(CDRRange(annotation.cdr1_start, annotation.cdr1_end),
                               CDRRange(annotation.cdr2_start, annotation.cdr2_end),
                               CDRRange(annotation.cdr3_start, size_t(-1)));
        }
        return CDRLabeling();
    }

    void AnnotatedJGeneCDRLabeler::Initialize() {
        std::string j_gene_annotation = search_params_.imgt_j_annotation;
        VERIFY_MSG(false, "Implement me!");
    }

    CDRLabeling AnnotatedJGeneCDRLabeler::ComputeLabeling(const germline_utils::ImmuneGene &immune_gene) {
        VERIFY_MSG(immune_gene.GeneType().Segment(), "Type of immune gene " << immune_gene.name() << " is not J");
        VERIFY_MSG(false, "Implement me!");
        return CDRLabeling();
    }

    CDRLabeling DeNovoImmuneGeneCDRLabeler::ComputeLabeling(const germline_utils::ImmuneGene &immune_gene) {
        VERIFY_MSG(immune_gene.GeneType() == gene_type_, "Type of immune gene (" << immune_gene.GeneType() <<
                                                         ") is not " << gene_type_);
        CDRRange cdr1 = cdr1_labeler_ptr_->ComputeLoopRange(immune_gene, CDRRange());
        CDRRange cdr2 = cdr2_labeler_ptr_->ComputeLoopRange(immune_gene, cdr1);
        CDRRange cdr3; // = cdr3_labeler_ptr_->ComputeLoopRange(immune_gene, cdr2);
        return CDRLabeling(cdr1, cdr2, cdr3);
    }

    SingleLoopLabelerPtr ImmuneGeneCDRHelper::GetCDR1Labeler(germline_utils::ImmuneGeneType gene_type,
                                                             const CDRLabelerConfig::CDRsParams &cdr_params) {
        if(gene_type.Segment() == germline_utils::SegmentType::JoinSegment)
            return SingleLoopLabelerPtr(new SingleLoopLabeler(gene_type));
        return SingleLoopLabelerPtr(new HCDR1Labeler(cdr_params.hcdr1_params));
    }

    SingleLoopLabelerPtr ImmuneGeneCDRHelper::GetCDR2Labeler(germline_utils::ImmuneGeneType gene_type,
                                                             const CDRLabelerConfig::CDRsParams &cdr_params) {
        if(gene_type.Segment() == germline_utils::SegmentType::JoinSegment)
            return SingleLoopLabelerPtr(new SingleLoopLabeler(gene_type));
        return SingleLoopLabelerPtr(new HCDR2Labeler(cdr_params.hcdr2_params));
    }

    SingleLoopLabelerPtr ImmuneGeneCDRHelper::GetCDR3Labeler(germline_utils::ImmuneGeneType gene_type,
                                                             const CDRLabelerConfig::CDRsParams &cdr_params) {
        if(gene_type.Segment() == germline_utils::SegmentType::JoinSegment)
            return SingleLoopLabelerPtr(new HCDR3JLabeler(cdr_params.hcdr3_params));
        return SingleLoopLabelerPtr(new HCDR3VLabeler(cdr_params.hcdr3_params));
    }
}
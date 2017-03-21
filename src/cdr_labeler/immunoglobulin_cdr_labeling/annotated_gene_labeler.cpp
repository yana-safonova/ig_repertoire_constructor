#include "annotated_gene_labeler.hpp"

#include <convert.hpp>

#include <boost/algorithm/string.hpp>

namespace cdr_labeler {
    using namespace annotation_utils;

    void AnnotatedVGeneCDRLabeler::Initialize() {
        auto v_annotation_params = search_params_.v_gene_annotation;
        std::string v_gene_annotation = (search_params_.domain_system == CDRLabelerConfig::CDRsParams::
        AnnotatedSearchParams::DomainSystem::IMGT_Domain) ?
                                        v_annotation_params.imgt_v_annotation : v_annotation_params.kabat_v_annotation;
        INFO("Reading V annotation from " << v_gene_annotation);
        VERIFY_MSG(path::check_existence(v_gene_annotation), "File with V gene annotation " << v_gene_annotation <<
                                                             " does not exist");
        std::ifstream ifhandler(v_gene_annotation);
        while(!ifhandler.eof()) {
            std::string tmp;
            std::getline(ifhandler, tmp);
            if(tmp == "")
                break;
            std::vector<std::string> splits;
            boost::split(splits, tmp, boost::is_any_of("\t "), boost::token_compress_on);
            AnnotatedVGeneCDRLabeler::VGeneAnnotation v_annotation(splits[v_annotation_params.v_gene_line_index],
                                                                   core::convert<std::string, size_t>(splits[v_annotation_params.cdr1_start_line_index]) - 1,
                                                                   core::convert<std::string, size_t>(splits[v_annotation_params.cdr1_end_line_index]) - 1,
                                                                   core::convert<std::string, size_t>(splits[v_annotation_params.cdr2_start_line_index]) - 1,
                                                                   core::convert<std::string, size_t>(splits[v_annotation_params.cdr2_end_line_index]) - 1,
                                                                   core::convert<std::string, size_t>(splits[v_annotation_params.fr3_end_index]));
            v_annotations_.push_back(v_annotation);
            v_name_index_map_[v_annotation.name] = v_annotations_.size() - 1;
        }
        //INFO(v_annotations_.size() << " V annotations were extracted from " << v_gene_annotation);
    }

    CDRLabeling AnnotatedVGeneCDRLabeler::ComputeLabeling(const germline_utils::ImmuneGene &immune_gene) {
        VERIFY_MSG(immune_gene.GeneType().Segment() == segment_type_, "Type of immune gene " << immune_gene.name() << " is not V");
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
        auto j_params = search_params_.j_gene_annotation;
        std::string j_gene_annotation = j_params.imgt_j_annotation;
        INFO("Reading J annotation from " << j_params.imgt_j_annotation);
        VERIFY_MSG(path::check_existence(j_gene_annotation), "File with J gene annotation " << j_gene_annotation <<
                                                             " does not exist");
        std::ifstream ifhandler(j_gene_annotation);
        while(!ifhandler.eof()) {
            std::string tmp;
            std::getline(ifhandler, tmp);
            if(tmp == "")
                break;
            std::vector<std::string> splits;
            boost::split(splits, tmp, boost::is_any_of("\t"));
            AnnotatedJGeneCDRLabeler::JGeneAnnotation j_annotation(splits[j_params.j_gene_line_index],
                                                                   core::convert<std::string, size_t>(splits[j_params.cdr3_end_index]) - 1);
            j_annotations_.push_back(j_annotation);
            j_name_index_map_[j_annotation.name] = j_annotations_.size() - 1;
        }
        //INFO(j_annotations_.size() << " J annotations were extracted from " << j_gene_annotation);
    }

    CDRLabeling AnnotatedJGeneCDRLabeler::ComputeLabeling(const germline_utils::ImmuneGene &immune_gene) {
        VERIFY_MSG(immune_gene.GeneType().Segment() == segment_type_, "Type of immune gene " << immune_gene.name() << " is not J");
        std::string gene_name = std::string(seqan::toCString(immune_gene.name()));
        if(j_name_index_map_.find(gene_name) != j_name_index_map_.end()) {
            AnnotatedJGeneCDRLabeler::JGeneAnnotation annotation = j_annotations_[j_name_index_map_.at(gene_name)];
            return CDRLabeling(CDRRange(), CDRRange(), CDRRange(size_t(-1), annotation.cdr3_end), annotation.orf);
        }
        return CDRLabeling();
    }
}
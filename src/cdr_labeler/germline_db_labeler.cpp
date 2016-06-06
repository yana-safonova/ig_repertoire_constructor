#include <verify.hpp>

#include "immunoglobulin_cdr_labeling/immune_gene_labeling_helper.hpp"
#include "germline_db_labeler.hpp"

namespace cdr_labeler {
    BaseImmuneGeneCDRLabelerPtr GermlineDbLabeler::GetImmuneGeneLabeler(germline_utils::ImmuneGeneType gene_type) {
        if(cdr_params_.cdr_search_algorithm == CDRLabelerConfig::CDRsParams::CDRSearchAlgorithm::AnnotatedSearch)
            return ImmuneGeneLabelingHelper::GetAnnotatedLabeler(gene_type.Segment(), cdr_params_);
        return BaseImmuneGeneCDRLabelerPtr(new DeNovoImmuneGeneCDRLabeler(gene_type, cdr_params_));
    }

    std::string cdr_search_algorithm_to_str(CDRLabelerConfig::CDRsParams::CDRSearchAlgorithm alg) {
        if(alg == CDRLabelerConfig::CDRsParams::CDRSearchAlgorithm::AnnotatedSearch)
            return "annotated";
        if(alg == CDRLabelerConfig::CDRsParams::CDRSearchAlgorithm::DeNovoSearch)
            return "de novo";
        return "";
    }

    DbCDRLabeling GermlineDbLabeler::ComputeLabeling() {
        DbCDRLabeling cdr_labeling(gene_db_);
        INFO("Algorithm of CDR computation: " << cdr_search_algorithm_to_str(cdr_params_.cdr_search_algorithm));
        for(auto it = gene_db_.cbegin(); it != gene_db_.cend(); it++) {
            auto specific_gene_db = gene_db_.GetDbByGeneType(*it);
            //INFO("Processing DB of type " << specific_gene_db.Chain());
            auto cdr_labeler = GetImmuneGeneLabeler(specific_gene_db.GeneType());
            for(auto it2 = specific_gene_db.cbegin(); it2 != specific_gene_db.cend(); it2++) {
                cdr_labeling.AddGeneLabeling(*it2, cdr_labeler->ComputeLabeling(*it2));
            }
        }
        INFO("# records from DB with empty CDR labelings: " << cdr_labeling.NumEmptyLabelings());
        return cdr_labeling;
    }
}
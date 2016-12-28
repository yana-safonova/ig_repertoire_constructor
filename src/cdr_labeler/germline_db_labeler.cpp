#include <verify.hpp>

#include "immunoglobulin_cdr_labeling/immune_gene_labeling_helper.hpp"
#include "germline_db_labeler.hpp"

namespace cdr_labeler {
    BaseImmuneGeneCDRLabelerPtr GermlineDbLabeler::GetImmuneGeneLabeler(germline_utils::ImmuneGeneType gene_type) {
        if(cdr_params_.cdr_search_algorithm == CDRLabelerConfig::CDRsParams::CDRSearchAlgorithm::AnnotatedCDRSearch)
            return ImmuneGeneLabelingHelper::GetAnnotatedLabeler(gene_type.Segment(), cdr_params_);
        return BaseImmuneGeneCDRLabelerPtr(new DeNovoImmuneGeneCDRLabeler(gene_type, cdr_params_));
    }

    std::string cdr_search_algorithm_to_str(CDRLabelerConfig::CDRsParams::CDRSearchAlgorithm alg) {
        if(alg == CDRLabelerConfig::CDRsParams::CDRSearchAlgorithm::AnnotatedCDRSearch)
            return "annotated";
        if(alg == CDRLabelerConfig::CDRsParams::CDRSearchAlgorithm::DeNovoCDRSearch)
            return "de novo";
        return "";
    }

    DbCDRLabeling GermlineDbLabeler::ComputeLabeling() {
        DbCDRLabeling cdr_labeling(gene_db_);
        INFO("Algorithm of CDR computation: " << cdr_search_algorithm_to_str(cdr_params_.cdr_search_algorithm));
        for(auto it = gene_db_.cbegin(); it != gene_db_.cend(); it++) {
            germline_utils::ImmuneGeneDatabase& specific_gene_db = gene_db_.GetDbByGeneType(*it);
            auto cdr_labeler = GetImmuneGeneLabeler(specific_gene_db.GeneType());
            for(size_t i = 0; i < specific_gene_db.size(); i++) {
                auto gene_labeling = cdr_labeler->ComputeLabeling(specific_gene_db[i]);
                specific_gene_db.GetImmuneGeneByIndex(i).SetORF(gene_labeling.orf);
                cdr_labeling.AddGeneLabeling(specific_gene_db[i], gene_labeling);
            }
        }
        INFO("# records from DB with empty CDR labelings: " << cdr_labeling.NumEmptyLabelings());
        return cdr_labeling;
    }
}
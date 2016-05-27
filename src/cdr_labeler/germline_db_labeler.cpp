#include <verify.hpp>

#include "immunoglobulin_cdr_labeling/immune_gene_labeler.hpp"
#include "germline_db_labeler.hpp"

namespace cdr_labeler {
    DbCDRLabeling GermlineDbLabeler::ComputeLabeling() {
        DbCDRLabeling cdr_labeling(gene_db_);
        for(auto it = gene_db_.cbegin(); it != gene_db_.cend(); it++) {
            auto specific_gene_db = gene_db_.GetDbByGeneType(*it);
            ImmuneGeneCDRLabeler gene_labeler(*it, cdr_params_);
            for(auto it2 = specific_gene_db.cbegin(); it2 != specific_gene_db.cend(); it2++) {
                std::cout << "Processing immune gene " << *it2 << std::endl;
                cdr_labeling.AddGeneLabeling(*it2, gene_labeler.ComputeLabeling(*it2));
            }
        }
        return cdr_labeling;
    }
}
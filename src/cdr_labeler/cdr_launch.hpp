#pragma once

#include "cdr_config.hpp"
#include "germline_db_labeling.hpp"

namespace cdr_labeler {
    class CDRLabelerLaunch {
        const CDRLabelerConfig &config_;

        germline_utils::CustomGeneDatabase GetDatabaseByCDRLabeling(const germline_utils::CustomGeneDatabase &gene_db,
                                                                    DbCDRLabeling cdr_labeling);

    public:
        CDRLabelerLaunch(const CDRLabelerConfig &config) :
                config_(config) { }

        void Launch();
    };
}
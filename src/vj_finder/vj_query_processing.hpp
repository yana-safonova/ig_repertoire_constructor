#pragma once

#include "vj_query_aligner.hpp"
#include "vj_alignment_structs.hpp"
#include "vj_finder_config.hpp"

#include <boost/optional.hpp>

namespace vj_finder {
    typedef boost::optional<VJHits> ProcessedVJHits;

    class VJQueryProcessor {
        const VJFinderConfig::AlgorithmParams &params_;
        const germline_utils::CustomGeneDatabase &v_db_;
        const germline_utils::CustomGeneDatabase &j_db_;

        ProcessedVJHits ComputeFilteringResults(VJHits vj_hits);

    public:
        VJQueryProcessor(const VJFinderConfig::AlgorithmParams &params,
                         const germline_utils::CustomGeneDatabase &v_db,
                         const germline_utils::CustomGeneDatabase &j_db) : params_(params),
                                                                            v_db_(v_db),
                                                                            j_db_(j_db) { }

        ProcessedVJHits Process(const core::Read &read);
    };
}
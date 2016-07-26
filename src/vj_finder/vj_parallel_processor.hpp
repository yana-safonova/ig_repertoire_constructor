#pragma once

#include "vj_alignment_info.hpp"
#include "vj_query_processing.hpp"

namespace vj_finder {
    class VJParallelProcessor {
        core::ReadArchive &read_archive_;
        const VJFinderConfig::AlgorithmParams &algorithm_params_;
        const germline_utils::CustomGeneDatabase &v_db_;
        const germline_utils::CustomGeneDatabase &j_db_;
        size_t num_threads_;

        // i-th element shows which thread processed i-th read
        std::vector<size_t> thread_id_per_read_;
        // i-th element stores Alignment info created by i-th thread
        std::vector<VJAlignmentInfo> info_per_thread;

        void Initialize();

        VJAlignmentInfo GatherAlignmentInfos();

    public:
        VJParallelProcessor(core::ReadArchive &read_archive,
                            const VJFinderConfig::AlgorithmParams &algorithm_params,
                            const germline_utils::CustomGeneDatabase &v_db,
                            const germline_utils::CustomGeneDatabase &j_db_,
                            size_t num_threads) : read_archive_(read_archive),
                                                  algorithm_params_(algorithm_params),
                                                  v_db_(v_db), j_db_(j_db_),
                                                  num_threads_(num_threads) {
            Initialize();
        }

        VJAlignmentInfo Process();
    };
}
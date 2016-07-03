#include "vj_parallel_processor.hpp"

namespace vj_finder {
    void VJParallelProcessor::Initialize() {
        for(size_t i = 0; i < read_archive_.size(); i++)
            thread_id_per_read_.push_back(size_t(-1));
        for(size_t i = 0; i < num_threads_; i++)
            info_per_thread.push_back(VJAlignmentInfo());
    }

    VJAlignmentInfo VJParallelProcessor::GatherAlignmentInfos() {
        VJAlignmentInfo consistent_alignment_info;
        for(size_t i = 0; i < read_archive_.size(); i++) {
            size_t thread_id = thread_id_per_read_[i];
            if(info_per_thread[thread_id].ReadIsFiltered(read_archive_[i]))
                consistent_alignment_info.UpdateFilteredReads(read_archive_[i]);
            else {
                size_t local_index = info_per_thread[thread_id].GetVJHitIndexByRead(read_archive_[i]);
                consistent_alignment_info.UpdateHits(info_per_thread[thread_id].GetVJHitsByIndex(local_index));
            }
        }
        return consistent_alignment_info;
    }

    VJAlignmentInfo VJParallelProcessor::Process() {
        omp_set_num_threads(int(num_threads_));
#pragma omp parallel for schedule(dynamic)
        for(size_t i = 0; i < read_archive_.size(); i++) {
            size_t thread_id = omp_get_thread_num();
            thread_id_per_read_[i] = thread_id;
            VJQueryProcessor vj_query_processor(algorithm_params_, v_db_, j_db_);
            auto vj_hits = vj_query_processor.Process(read_archive_[i]);
            if(!vj_hits)
                info_per_thread[thread_id].UpdateFilteredReads(read_archive_[i]);
            else
                info_per_thread[thread_id].UpdateHits(*vj_hits);
        }
        return GatherAlignmentInfos();
    }
}
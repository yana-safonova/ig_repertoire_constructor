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
            if(info_per_thread[thread_id].ReadIsFiltered(read_archive_[i])) {
                size_t local_index = info_per_thread[thread_id].GetFilteringInfoIndexByRead(read_archive_[i]);
                consistent_alignment_info.UpdateFilteringInfo(
                        info_per_thread[thread_id].GetFilteringInfoByIndex(local_index));
            }
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
            TRACE("Processing read: " << read_archive_[i].name);
            size_t thread_id = omp_get_thread_num();
            thread_id_per_read_[i] = thread_id;
            VJQueryProcessor vj_query_processor(algorithm_params_, read_archive_, v_db_, j_db_);
            auto processed_read = vj_query_processor.Process(read_archive_[i]);
            if(processed_read.ReadToBeFiltered()) {
//                std::cout << "bad: " << processed_read.filtering_info.filtering_reason << std::endl;
                info_per_thread[thread_id].UpdateFilteringInfo(processed_read.filtering_info);
            }
            else {
//                std::cout << "good" << std::endl;
                info_per_thread[thread_id].UpdateHits(processed_read.vj_hits);
            }
        }
        auto total_alignment_info = GatherAlignmentInfos();
        size_t num_aligned_reads = total_alignment_info.NumVJHits();
        for(auto it = total_alignment_info.chain_type_cbegin(); it != total_alignment_info.chain_type_cend(); it++) {
            float perc = float(it->second) / float(num_aligned_reads) * 100;
            INFO(perc << "% of aligned reads have isotype " << it->first);
        }
        return total_alignment_info;
    }
}

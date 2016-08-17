#pragma once

#include <unordered_set>
#include "vj_finder_config.hpp"
#include "vj_alignment_structs.hpp"
#include "vj_hits_filter.hpp"

namespace vj_finder {
    class VJAlignmentInfo {
        std::vector<VJHits> alignment_records_;
        std::vector<VJFilteringInfo> filtering_infos_;

        std::map<size_t, size_t> read_id_hit_index_map_;
        std::map<size_t, size_t> read_id_filtering_info_map_;

        std::unordered_map<germline_utils::ChainType, size_t, germline_utils::ChainTypeHasher> chain_type_abundance_;

        void UpdateChainTypeMap(const VJHits &vj_hits);

    public:
        VJAlignmentInfo() { }

        void UpdateHits(VJHits vj_hits);

        void UpdateFilteringInfo(VJFilteringInfo filtering_info);

        void Update(VJAlignmentInfo vj_alignment_info);

        size_t NumFilteredReads() const { return filtering_infos_.size(); }

        size_t NumVJHits() const { return alignment_records_.size(); }


        const VJHits& GetVJHitsByIndex(size_t index) const;

        const VJFilteringInfo GetFilteringInfoByIndex(size_t index) const;


        const VJHits& GetVJHitsByRead(const core::Read &read) const;

        VJFilteringInfo GetFilteringInfoByRead(const core::Read &read) const;

        bool ReadIsFiltered(const core::Read &read) const {
            return read_id_filtering_info_map_.find(read.id) != read_id_filtering_info_map_.end();
        }

        size_t GetVJHitIndexByRead(const core::Read &read) const {
            VERIFY_MSG(read_id_hit_index_map_.find(read.id) != read_id_hit_index_map_.end(),
                       "Info does contain record for aligned read " << read.name);
            return read_id_hit_index_map_.at(read.id);
        }

        size_t GetFilteringInfoIndexByRead(const core::Read &read) const {
            VERIFY_MSG(read_id_filtering_info_map_.find(read.id) != read_id_filtering_info_map_.end(),
                       "Info does contain record for filtered read " << read.name);
            return read_id_filtering_info_map_.at(read.id);
        }

        const core::Read& GetFilteredReadByIndex(size_t filtering_info_index) const {
            VERIFY_MSG(filtering_info_index < filtering_infos_.size(),
                       "Index exceeds number of filtering records");
            return *(filtering_infos_[filtering_info_index].read);
        }

        typedef std::unordered_map<germline_utils::ChainType, size_t,
                germline_utils::ChainTypeHasher>::const_iterator ChainTypeAbundanceConstIter;

        ChainTypeAbundanceConstIter chain_type_cbegin() const { return chain_type_abundance_.cbegin(); }

        ChainTypeAbundanceConstIter chain_type_cend() const { return chain_type_abundance_.cend(); }
    };

    class VJAlignmentOutput {
        const VJFinderConfig::IOParams::OutputParams &output_params_;
        const VJAlignmentInfo &alignment_info_;

    public:
        VJAlignmentOutput(const VJFinderConfig::IOParams::OutputParams &output_params,
                          const VJAlignmentInfo &alignment_info) :
                output_params_(output_params), alignment_info_(alignment_info) { }

        void OutputFilteredReads() const;

        void OutputAlignmentInfo() const;

        void OutputCleanedReads() const;

        void OutputVAlignments() const;

        void OutputFilteringInfo() const;
    };
}
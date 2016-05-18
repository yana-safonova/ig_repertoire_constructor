#include <verify.hpp>
#include "vj_alignment_info.hpp"

namespace vj_finder {
    void VJAlignmentInfo::Update(VJAlignmentInfo vj_alignment_info) {
        for(size_t i = 0; i < vj_alignment_info.NumVJHits(); i++)
            UpdateHits(vj_alignment_info.GetVJHitsByIndex(i));
        for(size_t i = 0; i < vj_alignment_info.NumFilteredReads(); i++)
            UpdateFilteredReads(vj_alignment_info.GetFilteredReadByIndex(i));
    }

    void VJAlignmentInfo::UpdateFilteredReads(const core::Read &read) {
        filtered_reads_.push_back(&read);
    }

    void VJAlignmentInfo::UpdateHits(VJHits vj_hits) {
        alignment_records_.push_back(std::move(vj_hits));
    }

    const core::Read& VJAlignmentInfo::GetFilteredReadByIndex(size_t index) const {
        VERIFY_MSG(index < NumFilteredReads(), "Index " << index << " exceeds number of filtered reads " <<
                NumFilteredReads());
        return *filtered_reads_[index];
    }

    const VJHits& VJAlignmentInfo::GetVJHitsByIndex(size_t index) const {
        VERIFY_MSG(index < NumVJHits(), "Index " << index << " exceeds number of VJ hits " <<
                                               NumVJHits());
        return alignment_records_[index];
    }

    void VJAlignmentOutput::OutputAlignmentInfo() const {
        std::ofstream out(output_params_.output_files.alignment_info_fname);
        for(size_t i = 0; i < alignment_info_.NumVJHits(); i++) {
            auto vj_hits = alignment_info_.GetVJHitsByIndex(i);
            for(size_t j = 0; j < output_params_.output_details.num_aligned_candidates; j++)
                out << vj_hits.Read().name << "\t" << vj_hits.GetVHitByIndex(j).Start() << "\t" <<
                        vj_hits.GetVHitByIndex(j).End() << "\t" <<
                        vj_hits.GetVHitByIndex(j).Score() << "\t" <<
                        vj_hits.GetVHitByIndex(j).ImmuneGene().name() << "\t" <<
                        vj_hits.GetJHitByIndex(j).Start() << "\t" <<
                        vj_hits.GetJHitByIndex(j).End() << "\t" << vj_hits.GetJHitByIndex(j).Score() << "\t" <<
                        vj_hits.GetJHitByIndex(j).ImmuneGene().name() << std::endl;
        }
        out.close();
    }

    void VJAlignmentOutput::OutputCleanedReads() const {
        std::ofstream out(output_params_.output_files.cleaned_reads_fname);
        for(size_t i = 0; i < alignment_info_.NumVJHits(); i++) {
            auto read = alignment_info_.GetVJHitsByIndex(i).Read();
            out << ">" << read.name << std::endl;
            out << read.seq << std::endl;
        }
        out.close();
    }

    void VJAlignmentOutput::OutputFilteredReads() const {
        std::ofstream out(output_params_.output_files.filtered_reads_fname);
        for(size_t i = 0; i < alignment_info_.NumFilteredReads(); i++) {
            auto read = alignment_info_.GetFilteredReadByIndex(i);
            out << ">" << read.name << std::endl;
            out << read.seq << std::endl;
        }
        out.close();
    }
}
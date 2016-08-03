#include <verify.hpp>
#include "vj_alignment_info.hpp"
#include "immune_gene_alignment_converter.hpp"

namespace vj_finder {
    void VJAlignmentInfo::Update(VJAlignmentInfo vj_alignment_info) {
        for(size_t i = 0; i < vj_alignment_info.NumVJHits(); i++)
            UpdateHits(vj_alignment_info.GetVJHitsByIndex(i));
        for(size_t i = 0; i < vj_alignment_info.NumFilteredReads(); i++)
            UpdateFilteredReads(vj_alignment_info.GetFilteredReadByIndex(i));
    }

    void VJAlignmentInfo::UpdateFilteredReads(const core::Read &read) {
        filtered_reads_.push_back(&read);
        filtered_read_ids_.insert(read.id);
    }

    void VJAlignmentInfo::UpdateChainTypeMap(const VJHits &vj_hits) {
        auto chain_type = vj_hits.GetVHitByIndex(0).ImmuneGene().Chain();
        if(chain_type_abundance_.find(chain_type) == chain_type_abundance_.end())
            chain_type_abundance_[chain_type] = 0;
        chain_type_abundance_[chain_type]++;
    }

    void VJAlignmentInfo::UpdateHits(VJHits vj_hits) {
        UpdateChainTypeMap(vj_hits);
        alignment_records_.push_back(std::move(vj_hits));
        read_id_hit_index_map_[vj_hits.Read().id] = alignment_records_.size() - 1;
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

    const VJHits& VJAlignmentInfo::GetVJHitsByRead(const core::Read &read) const {
        VERIFY_MSG(read_id_hit_index_map_.find(read.id) != read_id_hit_index_map_.end(),
                   "Alignment info does not contain read " << read.name);
        return alignment_records_[read_id_hit_index_map_.at(read.id)];
    }

    void VJAlignmentOutput::OutputAlignmentInfo() const {
        std::ofstream out(output_params_.output_files.alignment_info_fname);
        for(size_t i = 0; i < alignment_info_.NumVJHits(); i++) {
            auto vj_hits = alignment_info_.GetVJHitsByIndex(i);
            for(size_t j = 0; j < output_params_.output_details.num_aligned_candidates; j++)
                out << vj_hits.Read().name << "\t" << vj_hits.GetVHitByIndex(j).Start() + 1 << "\t" <<
                        vj_hits.GetVHitByIndex(j).End() << "\t" <<
                        vj_hits.GetVHitByIndex(j).Score() << "\t" <<
                        vj_hits.GetVHitByIndex(j).ImmuneGene().name() << "\t" <<
                        vj_hits.GetJHitByIndex(j).Start() + 1 << "\t" <<
                        vj_hits.GetJHitByIndex(j).End() << "\t" << vj_hits.GetJHitByIndex(j).Score() << "\t" <<
                        vj_hits.GetJHitByIndex(j).ImmuneGene().name() << std::endl;
        }
        out.close();
        INFO("Alignment info was written to " << output_params_.output_files.alignment_info_fname);
    }

    void VJAlignmentOutput::OutputCleanedReads() const {
        std::ofstream out(output_params_.output_files.cleaned_reads_fname);
        for(size_t i = 0; i < alignment_info_.NumVJHits(); i++) {
            auto read = alignment_info_.GetVJHitsByIndex(i).Read();
            out << ">" << read.name << std::endl;
            out << read.seq << std::endl;
        }
        out.close();
        INFO("Cleaned reads were written to " << output_params_.output_files.cleaned_reads_fname);
    }

    void VJAlignmentOutput::OutputFilteredReads() const {
        std::ofstream out(output_params_.output_files.filtered_reads_fname);
        for(size_t i = 0; i < alignment_info_.NumFilteredReads(); i++) {
            auto read = alignment_info_.GetFilteredReadByIndex(i);
            out << ">" << read.name << std::endl;
            out << read.seq << std::endl;
        }
        out.close();
        INFO("Cleaned reads were written to " << output_params_.output_files.filtered_reads_fname);
    }

    void VJAlignmentOutput::OutputVAlignments() const {
        ImmuneGeneAlignmentConverter alignment_converter;
        std::ofstream out(output_params_.output_files.valignments_filename);
        for(size_t i = 0; i < alignment_info_.NumVJHits(); i++) {
            auto vj_hits = alignment_info_.GetVJHitsByIndex(i);
            auto v_hit = vj_hits.GetVHitByIndex(0);
            auto v_alignment = alignment_converter.ConvertToAlignment(v_hit.ImmuneGene(), vj_hits.Read(),
                                                                      v_hit.BlockAlignment());
            auto subject_row = seqan::row(v_alignment.Alignment(), 0);
            auto query_row = seqan::row(v_alignment.Alignment(), 1);
            out << ">INDEX:" << i + 1 << "|READ:" << vj_hits.Read().name << "|START_POS:" <<
                    v_alignment.StartSubjectPosition() << "|END_POS:" <<
                    v_alignment.EndSubjectPosition() << std::endl;
            out << query_row << std::endl;
            out << ">INDEX:" << i + 1 << "|GENE:" << v_alignment.subject().name() <<
            "|START_POS:" << v_alignment.StartQueryPosition() << "|END_POS:" <<
                    v_alignment.EndQueryPosition() << "|CHAIN_TYPE:" <<
                    v_alignment.subject().Chain() << std::endl;
            out << subject_row << std::endl;
        }
        out.close();
        INFO("V alignments were written to " << output_params_.output_files.valignments_filename);
    }
}
#include "shm_calculator.hpp"

namespace annotation_utils {
    // todo: remove temporary stub
    char get_aa_by_pos(const AAString &aa, size_t nucl_pos, size_t orf) {
        if(nucl_pos < orf)
            return '-';
        size_t aa_pos = (nucl_pos - orf) / 3;
        if (aa_pos < seqan::length(aa))
            return aa[aa_pos];
        else if (aa_pos == seqan::length(aa))
            return '-';
        VERIFY_MSG(false, "Unknown position " << aa_pos << " of amino acid sequence " << aa);
        return '-';
    }

    GeneSegmentSHMs NaiveSHMCalculator::ComputeSHMs(const alignment_utils::ImmuneGeneReadAlignment& alignment,
                                                    const AminoAcidAnnotation<core::Read>& aa_annotation,
                                                    const CDRLabeling &) {
        GeneSegmentSHMs shms(alignment.query(), alignment.subject());
        auto gene_row = seqan::row(alignment.Alignment(), 0);
        auto read_row = seqan::row(alignment.Alignment(), 1);
        for(size_t i = alignment.RealStartAlignmentPos(); i <= alignment.RealEndAlignmentPos(); i++) {
            if(gene_row[i] != read_row[i]) {
                size_t real_read_pos = seqan::toSourcePosition(read_row, i);
                size_t real_gene_pos = seqan::toSourcePosition(gene_row, i);
                SHM shm(alignment.subject().Segment(), real_gene_pos, real_read_pos, gene_row[i], read_row[i],
                        get_aa_by_pos(alignment.subject().aa_seq(), real_gene_pos, alignment.subject().ORF()),
                        aa_annotation.GetAminoAcidByPos(real_read_pos));
                shms.AddSHM(shm);
            }
        }
        return shms;
    }

    //--------------------------------------------------------------------

    void StartEndFilteringSHMCalculator::ComputeStartMeaningPositions(const GeneSegmentSHMs &all_shms,
            size_t, size_t) {
        if(all_shms.size() == 0) {
            //std::cout << "SHMs are empty" << std::endl;
            return;
        }
        first_meaning_read_pos_ = all_shms[0].read_nucl_pos;
        first_meaning_gene_pos_ = all_shms[0].gene_nucl_pos;
        for(size_t i = 1; i < all_shms.size(); i++) {
            size_t gene_diff = all_shms[i].gene_nucl_pos - first_meaning_gene_pos_;
            size_t read_diff = all_shms[i].read_nucl_pos - first_meaning_read_pos_;
            if(all_shms[i].gene_nucl_pos <= max_skipped_start_ or gene_diff <= 1 or read_diff <= 1) {
                first_meaning_gene_pos_ = all_shms[i].gene_nucl_pos;
                first_meaning_read_pos_ = all_shms[i].read_nucl_pos;
            }
            else
                break;
        }
        //std::cout << "First: " << first_meaning_read_pos_ << ", " << first_meaning_gene_pos_ << std::endl;
        first_meaning_read_pos_ = std::max(max_skipped_start_, first_meaning_read_pos_);
        first_meaning_gene_pos_ = std::max(max_skipped_start_, first_meaning_gene_pos_);
    }

    void StartEndFilteringSHMCalculator::ComputeEndMeaningPositions(const GeneSegmentSHMs &all_shms,
                                                                    size_t end_read_pos, size_t end_gene_pos) {
        if(all_shms.size() == 0) {
            //std::cout << "SHMs are empty" << std::endl;
            return;
        }
        last_meaning_read_pos_ = all_shms[all_shms.size() - 1].read_nucl_pos;
        last_meaning_gene_pos_ = all_shms[all_shms.size() - 1].gene_nucl_pos;
        for(size_t i = 1; i < all_shms.size(); i++) {
            size_t gene_diff = last_meaning_gene_pos_ - all_shms[all_shms.size() - i - 1].gene_nucl_pos;
            size_t read_diff = last_meaning_read_pos_ - all_shms[all_shms.size() - i - 1].read_nucl_pos;
            size_t length_from_gene_end = all_shms.ImmuneGene().length() - all_shms[all_shms.size() - i - 1].gene_nucl_pos;
            if(length_from_gene_end <= max_skipped_end_ or gene_diff <= 1 or read_diff <= 1) {
                last_meaning_gene_pos_ = all_shms[all_shms.size() - i - 1].gene_nucl_pos;
                last_meaning_read_pos_ = all_shms[all_shms.size() - i - 1].read_nucl_pos;
            }
            else
                break;
        }
        last_meaning_read_pos_ = std::min(end_read_pos - max_skipped_end_,
                                          last_meaning_read_pos_);
        last_meaning_gene_pos_ = std::min(end_gene_pos - max_skipped_end_,
                                          last_meaning_gene_pos_);
    }

    void StartEndFilteringSHMCalculator::ComputeMeaningPositions(const GeneSegmentSHMs& all_shms,
                                                                 const alignment_utils::ImmuneGeneReadAlignment& alignment) {
        ComputeStartMeaningPositions(all_shms, alignment.StartQueryPosition(), alignment.StartSubjectPosition());
        ComputeEndMeaningPositions(all_shms, alignment.EndQueryPosition(), alignment.EndSubjectPosition());
    }

    GeneSegmentSHMs StartEndFilteringSHMCalculator::ComputeSHMs(
            const alignment_utils::ImmuneGeneReadAlignment &alignment,
            const AminoAcidAnnotation<core::Read> &aa_annotation,
            const CDRLabeling &cdr_labeling) {
        GeneSegmentSHMs all_shms = NaiveSHMCalculator().ComputeSHMs(alignment, aa_annotation, cdr_labeling);
        //std::cout << alignment.subject().Segment() << std::endl;
        //std::cout << all_shms << std::endl;
        ComputeMeaningPositions(all_shms, alignment);
        GeneSegmentSHMs filtered_shms(alignment.query(), alignment.subject());
        //std::cout << first_meaning_read_pos_ << " - " << first_meaning_gene_pos_ << std::endl;
        //std::cout << last_meaning_read_pos_ << " - " << last_meaning_gene_pos_ << std::endl;
        for(auto it = all_shms.cbegin(); it != all_shms.cend(); it++) {
            if(it->read_nucl_pos >= first_meaning_read_pos_ and it->gene_nucl_pos >= first_meaning_gene_pos_ and
                    it->read_nucl_pos <= last_meaning_read_pos_ and it->gene_nucl_pos <= last_meaning_gene_pos_) {
                filtered_shms.AddSHM(*it);
            }
        }
        return filtered_shms;
    }

    //--------------------------------------------------------------------

    void CDRFilteringSHMCalculator::ComputeStartMeaningPositions(const alignment_utils::ImmuneGeneReadAlignment &alignment,
                                                                 const GeneSegmentSHMs &all_shms,
                                                                 const CDRLabeling &cdr_labeling) {
        if (all_shms.size() == 0) {
            return;
        }
        if (all_shms.SegmentType() == germline_utils::SegmentType::VariableSegment or !cdr_labeling.cdr3.Valid()) {
            first_meaning_read_pos_ = all_shms[0].read_nucl_pos;
            first_meaning_gene_pos_ = all_shms[0].gene_nucl_pos;
            return;
        }
        VERIFY_MSG(all_shms.SegmentType() == germline_utils::SegmentType::JoinSegment,
                   "Segment " << all_shms.SegmentType() << " is not variable or diversity");
        VERIFY_MSG(cdr_labeling.cdr3.Valid(), "CDR3 is not defined");
        first_meaning_read_pos_ = cdr_labeling.cdr3.end_pos + 1;
        first_meaning_gene_pos_ = alignment.SubjectPositionByQueryPosition(first_meaning_read_pos_);
    }

    void CDRFilteringSHMCalculator::ComputeEndMeaningPositions(const alignment_utils::ImmuneGeneReadAlignment &alignment,
                                                               const GeneSegmentSHMs &all_shms,
                                                               const CDRLabeling &cdr_labeling) {
        if(all_shms.size() == 0) {
            return;
        }
        if(all_shms.SegmentType() == germline_utils::SegmentType::JoinSegment or !cdr_labeling.cdr3.Valid()) {
            last_meaning_gene_pos_ = all_shms[all_shms.size() - 1].gene_nucl_pos;
            last_meaning_read_pos_ = all_shms[all_shms.size() - 1].read_nucl_pos;
            return;
        }
        VERIFY_MSG(all_shms.SegmentType() == germline_utils::SegmentType::VariableSegment,
                   "Segment " << all_shms.SegmentType() << " is not variable or diversity");
        VERIFY_MSG(cdr_labeling.cdr3.Valid(), "CDR3 is not defined");
        last_meaning_read_pos_ = cdr_labeling.cdr3.start_pos - 1;
        last_meaning_gene_pos_ = alignment.SubjectPositionByQueryPosition(last_meaning_read_pos_);
    }

    void CDRFilteringSHMCalculator::ComputeMeaningPositions(const alignment_utils::ImmuneGeneReadAlignment &alignment,
                                                            const GeneSegmentSHMs &all_shms,
                                                            const CDRLabeling &cdr_labeling) {
        if(all_shms.size() == 0)
            return;
        ComputeStartMeaningPositions(alignment, all_shms, cdr_labeling);
        ComputeEndMeaningPositions(alignment, all_shms, cdr_labeling);
    }

    GeneSegmentSHMs CDRFilteringSHMCalculator::ComputeSHMs(const alignment_utils::ImmuneGeneReadAlignment &alignment,
                                                           const AminoAcidAnnotation<core::Read> &aa_annotation,
                                                           const CDRLabeling &cdr_labeling) {
        GeneSegmentSHMs all_shms = NaiveSHMCalculator().ComputeSHMs(alignment, aa_annotation, cdr_labeling);
        ComputeMeaningPositions(alignment, all_shms, cdr_labeling);
        GeneSegmentSHMs filtered_shms(alignment.query(), alignment.subject());
        for(auto it = all_shms.cbegin(); it != all_shms.cend(); it++) {
            if(it->read_nucl_pos >= first_meaning_read_pos_ and it->gene_nucl_pos >= first_meaning_gene_pos_ and
               it->read_nucl_pos <= last_meaning_read_pos_ and it->gene_nucl_pos <= last_meaning_gene_pos_) {
                filtered_shms.AddSHM(*it);
            }
        }
        return filtered_shms;
    }
}
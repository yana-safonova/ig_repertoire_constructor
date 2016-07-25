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
                                                    const AminoAcidAnnotation<core::Read>& aa_annotation) {
        GeneSegmentSHMs shms(alignment.query(), alignment.subject());
        auto gene_row = seqan::row(alignment.Alignment(), 0);
        auto read_row = seqan::row(alignment.Alignment(), 1);
        for(size_t i = alignment.RealStartAlignmentPos(); i <= alignment.RealEndAlignmentPos(); i++) {
            if(gene_row[i] != read_row[i]) {
                size_t real_read_pos = seqan::toSourcePosition(read_row, i);
                size_t real_gene_pos = seqan::toSourcePosition(gene_row, i);
                SHM shm(real_gene_pos, real_read_pos, gene_row[i], read_row[i],
                        get_aa_by_pos(alignment.subject().aa_seq(), real_gene_pos, alignment.subject().ORF()),
                        aa_annotation.GetAminoAcidByPos(real_read_pos));
                shms.AddSHM(shm);
            }
        }
        return shms;
    }

    //--------------------------------------------------------------------

    void StartEndFilteringSHMCalculator::ComputeStartMeaningPositions(const GeneSegmentSHMs &all_shms,
            size_t start_read_pos, size_t start_gene_pos) {
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
            const AminoAcidAnnotation<core::Read> &aa_annotation) {
        GeneSegmentSHMs all_shms = NaiveSHMCalculator().ComputeSHMs(alignment, aa_annotation);
        //std::cout << alignment.subject().Segment() << std::endl;
        //std::cout << all_shms << std::endl;
        ComputeMeaningPositions(all_shms, alignment);
        GeneSegmentSHMs filtered_shms(alignment.query(), alignment.subject());
        //std::cout << first_meaning_read_pos_ << " - " << first_meaning_gene_pos_ << std::endl;
        //std::cout << last_meaning_read_pos_ << " - " << last_meaning_gene_pos_ << std::endl;
        for(auto it = all_shms.cbegin(); it != all_shms.cend(); it++) {
            if(it->read_nucl_pos >= first_meaning_read_pos_ and it->gene_nucl_pos >= first_meaning_gene_pos_ and
                    it->read_nucl_pos <= last_meaning_read_pos_ and it->gene_nucl_pos <= last_meaning_gene_pos_)
                filtered_shms.AddSHM(*it);
        }
        return filtered_shms;
    }
}
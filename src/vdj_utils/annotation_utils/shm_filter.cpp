#include "shm_filter.hpp"

namespace annotation_utils {
    void PositionalSHMFilter::ComputeStartMeaningPositions() {
        if(shms_.size() == 0)
            return;
        first_meaning_read_pos_ = shms_[0].read_nucl_pos;
        first_meaning_gene_pos_ = shms_[0].gene_nucl_pos;
        for(size_t i = 1; i < shms_.size(); i++) {
            size_t gene_diff = shms_[i].gene_nucl_pos - first_meaning_gene_pos_;
            size_t read_diff = shms_[i].read_nucl_pos - first_meaning_read_pos_;
            if(shms_[i].gene_nucl_pos <= max_num_skipped_start_ or gene_diff <= 1 or read_diff <= 1) {
                first_meaning_gene_pos_ = shms_[i].gene_nucl_pos;
                first_meaning_read_pos_ = shms_[i].read_nucl_pos;
            }
            else
                break;
        }
        first_meaning_read_pos_ = std::max(max_num_skipped_start_, first_meaning_read_pos_);
        first_meaning_gene_pos_ = std::max(max_num_skipped_start_, first_meaning_gene_pos_);
    }

    void PositionalSHMFilter::ComputeEndMeaningPositions() {
        if(shms_.size() == 0)
            return;
        last_meaning_read_pos_ = shms_[shms_.size() - 1].read_nucl_pos;
        last_meaning_gene_pos_ = shms_[shms_.size() - 1].gene_nucl_pos;
        for(size_t i = 1; i < shms_.size(); i++) {
            size_t gene_diff = last_meaning_gene_pos_ - shms_[shms_.size() - i - 1].gene_nucl_pos;
            size_t read_diff = last_meaning_read_pos_ - shms_[shms_.size() - i - 1].read_nucl_pos;
            size_t length_from_gene_end = shms_.ImmuneGene().length() - shms_[shms_.size() - i - 1].gene_nucl_pos;
            if(length_from_gene_end <= max_num_skipped_end_ or gene_diff <= 1 or read_diff <= 1) {
                last_meaning_gene_pos_ = shms_[shms_.size() - i - 1].gene_nucl_pos;
                last_meaning_read_pos_ = shms_[shms_.size() - i - 1].read_nucl_pos;
            }
            else
                break;
        }
        last_meaning_read_pos_ = std::max(max_num_skipped_end_, last_meaning_read_pos_);
        last_meaning_gene_pos_ = std::max(max_num_skipped_end_, last_meaning_gene_pos_);
    }

    void PositionalSHMFilter::ComputeMeaningPositions() {
        ComputeStartMeaningPositions();
        ComputeEndMeaningPositions();
        //std::cout << shms_ << std::endl;
        //std::cout << first_meaning_read_pos_ << " - " << first_meaning_gene_pos_ << std::endl;
        //std::cout << last_meaning_read_pos_ << " - " << last_meaning_gene_pos_ << std::endl;
    }

    bool PositionalSHMFilter::FilterSHM(SHM shm) const {
        return !(shm.read_nucl_pos >= first_meaning_read_pos_ and shm.gene_nucl_pos >= first_meaning_gene_pos_ and
                shm.read_nucl_pos <= last_meaning_read_pos_ and shm.gene_nucl_pos <= last_meaning_gene_pos_);
    }
}
//
// Created by Andrew Bzikadze on 5/18/16.
//

#pragma once

#include "verify.hpp"

#include <string>
#include <vector>

#include "convert.hpp"
#include "annotation_utils/annotated_clone.hpp"

namespace shm_kmer_matrix_estimator {

// using DnaGapped = seqan::ModifiedAlphabet<seqan::Dna5, seqan::ModExpand<'-'>>;
// using DnaGappedString = seqan::String<DnaGapped>;
// using DnaGappedAlignment = seqan::Align<DnaGappedString, seqan::ArrayGaps>;

// We place germline and read as pure c++ strings here for now.
// In AntEvolo we need to calculate weight not only of germline->read edges but also of evolutionary tree's edges.
// Due to that we call the class EvolutionaryEdgeAlignment instead of GermlineReadAlignemnt
class EvolutionaryEdgeAlignment {
private:
    // first component is parent. It could be gene.
    // second component is son.
    std::string parent_;
    std::string son_;

    // id of a gene each string belongs to.
    std::string gene_id_;

    bool has_stop_codon_;
    bool in_frame_;
    bool productive_;

    size_t cdr1_start_;
    size_t cdr1_end_;
    size_t cdr2_start_;
    size_t cdr2_end_;


    bool cropped_;
    bool checked_;
    bool check_result_;

public:
    EvolutionaryEdgeAlignment(const std::string &parent,
                              const std::string &son,
                              const std::string &gene_id,
                              const bool has_stop_codon,
                              const bool in_frame,
                              const bool productive,
                              const size_t cdr1_start,
                              const size_t cdr1_end,
                              const size_t cdr2_start,
                              const size_t cdr2_end) :
        parent_(parent),
        son_(son),
        gene_id_(gene_id),
        has_stop_codon_(has_stop_codon),
        in_frame_(in_frame),
        productive_(productive),
        cdr1_start_(cdr1_start),
        cdr1_end_(cdr1_end),
        cdr2_start_(cdr2_start),
        cdr2_end_(cdr2_end),
        cropped_(false),
        checked_(false),
        check_result_(false)
    {
        VERIFY_MSG(parent.size() == son.size(),
                   "Parent and son lengths are not equal\n" + parent + "\n" + son + "\n" + gene_id);
        VERIFY_MSG(cdr1_start <= cdr1_end and cdr1_end <= cdr2_start and cdr2_start <= cdr2_end,
                   std::string("Incorrect cdr initial values: ") +
                   std::string("cdr1_start = ") + std::to_string(cdr1_start) +
                   std::string("cdr1_end   = ") + std::to_string(cdr1_end) +
                   std::string("cdr2_start = ") + std::to_string(cdr2_start) +
                   std::string("cdr2_end   = ") + std::to_string(cdr2_end));
    }

    EvolutionaryEdgeAlignment(const annotation_utils::AnnotatedClone& clone) :
        parent_(core::seqan_string_to_string(seqan::row(clone.VAlignment().Alignment(), 0))),
        son_   (core::seqan_string_to_string(seqan::row(clone.VAlignment().Alignment(), 1))),
        gene_id_(core::seqan_string_to_string(clone.VAlignment().subject().name())),
        has_stop_codon_(clone.HasStopCodon()),
        in_frame_(clone.InFrame()),
        productive_(clone.Productive()),
        cdr1_start_(clone.GetRangeByRegion(annotation_utils::StructuralRegion::CDR1).start_pos),
        cdr1_end_  (clone.GetRangeByRegion(annotation_utils::StructuralRegion::CDR1).end_pos),
        cdr2_start_(clone.GetRangeByRegion(annotation_utils::StructuralRegion::CDR2).start_pos),
        cdr2_end_  (clone.GetRangeByRegion(annotation_utils::StructuralRegion::CDR2).end_pos)
    {
        VERIFY_MSG(parent_.size() == son_.size(),
                   "Parent and son lengths are not equal\n" + parent_ + "\n" + son_ + "\n" + gene_id_);
        VERIFY_MSG(cdr1_start_ <= cdr1_end_ and cdr1_end_ <= cdr2_start_ and cdr2_start_ <= cdr2_end_,
                   "Incorrect cdr initial values");
    }

    const std::string &parent() const { return parent_; }
    const std::string &son() const { return son_; }
    const std::string &gene_id() const { return gene_id_; }

    size_t cdr1_start() const { return cdr1_start_; }
    size_t cdr2_start() const { return cdr2_start_; }

    size_t cdr1_end() const { return cdr1_end_; }
    size_t cdr2_end() const { return cdr2_end_; }

    void set_parent(const std::string &parent) { parent_ = parent; }
    void set_son(const std::string &son) { son_ = son; }

    void set_parent(const std::string::const_iterator &begin,
                    const std::string::const_iterator &end) {
        parent_.assign<std::string::const_iterator>(begin, end);
    }

    void set_son(const std::string::const_iterator &begin,
                 const std::string::const_iterator &end) {
        son_.assign<std::string::const_iterator>(begin, end);
    }

    size_t size() const {
        VERIFY_MSG(parent_.size() == son_.size(),
                   "Parent and son lengths are not equal\n" + parent_ + "\n" + son_ + "\n" + gene_id_);
        return parent_.size();
    }

    void substract_cdr_positions(size_t x) {
        VERIFY_MSG(cdr1_start_ >= x, "Incorrect substract value");
        cdr1_start_ -= x, cdr1_end_ -= x;
        cdr2_start_ -= x, cdr2_end_ -= x;
    }

    bool HasStopCodon() const { return has_stop_codon_; }
    bool InFrame() const { return in_frame_; }
    bool Productive() const { return productive_; }

    bool IsCropped() const { return cropped_; }
    bool IsChecked() const { return checked_; }
    bool CheckIsOk() const { return check_result_; }
    void SetCropped() { cropped_ = true; }
    void SetChecked() { checked_ = true; }
    void SetCheckResult(bool check_result) { check_result_ = check_result; }
};
using VectorEvolutionaryEdgeAlignments = std::vector<EvolutionaryEdgeAlignment>;

} // End namespace shm_kmer_matrix_estimator

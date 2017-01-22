//
// Created by Andrew Bzikadze on 5/18/16.
//

#pragma once

#include "verify.hpp"

#include <string>
#include <vector>

namespace ns_gene_alignment {
// using DnaGapped = seqan::ModifiedAlphabet<seqan::Dna5, seqan::ModExpand<'-'>>;
// using DnaGappedString = seqan::String<DnaGapped>;
// using DnaGappedAlignment = seqan::Align<DnaGappedString, seqan::ArrayGaps>;

// We place germline and read as pure c++ strings here for now.
// In AntEvolo we need to calculate weight not only of germline->read edges but also of evolutionary tree's edges.
// Due to that we call the class EvolutionaryEdgeAlignment instead of GermlineReadAlignemnt
class EvolutionaryEdgeAlignment {
private:
private:
    // first component is parent. It could be gene.
    // second component is son.
    std::string parent_;
    std::string son_;

    // id of a gene each string belongs to.
    std::string gene_id_;
    size_t cdr1_start_;
    size_t cdr1_end_;
    size_t cdr2_start_;
    size_t cdr2_end_;

public:
    EvolutionaryEdgeAlignment(const std::string &parent,
                              const std::string &son,
                              const std::string &gene_id,
                              const size_t cdr1_start = size_t(),
                              const size_t cdr1_end = size_t(),
                              const size_t cdr2_start = size_t(),
                              const size_t cdr2_end = size_t()) :
        parent_(parent),
        son_(son),
        gene_id_(gene_id),
        cdr1_start_(cdr1_start),
        cdr1_end_(cdr1_end),
        cdr2_start_(cdr2_start),
        cdr2_end_(cdr2_end)
    {
        VERIFY(parent.size() == son.size());
        VERIFY(cdr1_start <= cdr1_end <= cdr2_start <= cdr2_end);
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
        VERIFY(parent_.size() == son_.size());
        return parent_.size();
    }

    void substract_cdr_positions(size_t x) {
        VERIFY(cdr1_start_ >= x);
        cdr1_start_ -= x, cdr1_end_ -= x;
        cdr2_start_ -= x, cdr2_end_ -= x;
    }
};
using VectorEvolutionaryEdgeAlignments = std::vector<EvolutionaryEdgeAlignment>;
}

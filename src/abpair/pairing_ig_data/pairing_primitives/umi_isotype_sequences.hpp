#pragma once

#include "ig_isotype.hpp"
#include "seqan/sequence.h"

struct IsotypeUmiSequence {
    IgIsotype isotype;
    std::string umi;
    size_t size;
    seqan::Dna5String sequence;

    // temporary stub
    std::string v_gene;
    std::string j_gene;

    IsotypeUmiSequence() :
            isotype(),
            umi(),
            size(),
            sequence() { }

    IsotypeUmiSequence(IgIsotype new_isotype,
                       std::string new_umi,
                       size_t new_size,
                       seqan::Dna5String new_sequence) :
            isotype(new_isotype),
            umi(new_umi),
            size(new_size),
            sequence(new_sequence) { }

};

//--------------------------------------------------------------

class IsotypeUmiSequences {
    IgIsotype isotype_;
    std::vector<IsotypeUmiSequence> sequences_;

    void CheckConsistency(const IsotypeUmiSequence &umi_sequence) const;

public:
    IsotypeUmiSequences() : isotype_() { }

    IsotypeUmiSequences(IgIsotype isotype) : isotype_(isotype) { }

    void Update(IsotypeUmiSequence umi_sequence);

    size_t size() const { return sequences_.size(); }

    typedef std::vector<IsotypeUmiSequence>::const_iterator umi_seq_citerator;

    umi_seq_citerator cbegin() const { return sequences_.cbegin(); }

    umi_seq_citerator cend() const { return sequences_.cend(); }

    IgIsotype Isotype() const { return isotype_; }
};

typedef std::shared_ptr<IsotypeUmiSequences> IsotypeUmiSequencesPtr;
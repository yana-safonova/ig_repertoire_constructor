#include "umi_isotype_sequences.hpp"

void IsotypeUmiSequences::CheckConsistency(const IsotypeUmiSequence &umi_sequence) const {
    assert(isotype_ == umi_sequence.isotype);
}

void IsotypeUmiSequences::Update(IsotypeUmiSequence umi_sequence) {
    CheckConsistency(umi_sequence);
    sequences_.push_back(umi_sequence);
}
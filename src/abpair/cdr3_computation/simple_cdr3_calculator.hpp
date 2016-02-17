#pragma once

#include <vector>

#include "../pairing_ig_data/raw_pairing_data_storage.hpp"
#include <seqan/sequence.h>

class SimpleCdr3Calculator {
    bool PositionIsCys(const seqan::Dna5String &seq, size_t pos);

    bool PositionIsPhe(const seqan::Dna5String &seq, size_t pos);

    bool PositionIsTrp(const seqan::Dna5String &seq, size_t pos);

    std::vector<size_t> ComputeAaPositions(const seqan::Dna5String &seq, std::string aa);

    bool Cdr3PositionsCanBeRefined(std::pair<size_t, size_t> new_pos, std::pair<size_t, size_t> old_pos);

public:
    std::pair<size_t, size_t> FindCdr3Positions(const seqan::Dna5String &seq);
};
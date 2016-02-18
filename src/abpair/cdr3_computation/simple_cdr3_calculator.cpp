#include <logger/logger.hpp>
#include <verify.hpp>
#include "simple_cdr3_calculator.hpp"
#include <seqan/seq_io.h>

size_t abs_diff(size_t a, size_t b) {
    if(a > b)
        return a - b;
    return b - a;
}

bool SimpleCdr3Calculator::PositionIsCys(const seqan::Dna5String &seq, size_t pos) {
    VERIFY_MSG(pos + 2 < seqan::length(seq), "Cys cannot be found at last and prelast positions");
    return seq[pos] == 'T' and seq[pos + 1] == 'G' and seq[pos + 2] == 'T'; //(seq[pos + 2] == 'T' or seq[pos + 2] == 'C');
}

bool SimpleCdr3Calculator::PositionIsPhe(const seqan::Dna5String &seq, size_t pos) {
    VERIFY_MSG(pos + 2 < seqan::length(seq), "Phe cannot be found at last and prelast positions");
    return seq[pos] == 'T' and seq[pos + 1] == 'T' and (seq[pos + 2] == 'T' or seq[pos + 2] == 'C');
}

bool SimpleCdr3Calculator::PositionIsTrp(const seqan::Dna5String &seq, size_t pos) {
    VERIFY_MSG(pos + 2 < seqan::length(seq), "Trp cannot be found at last and prelast positions");
    return seq[pos] == 'T' and seq[pos + 1] == 'G' and seq[pos + 2] == 'G';
}

std::vector<size_t> SimpleCdr3Calculator::ComputeAaPositions(const seqan::Dna5String &seq, std::string aa) {
    std::vector<size_t> aa_positions;
    for(size_t i = 0; i < seqan::length(seq) - 2; i++) {
        bool pos_is_good = false;
        if(aa == "Phe")
            pos_is_good = PositionIsPhe(seq, i);
        else if(aa == "Cys")
            pos_is_good = PositionIsCys(seq, i);
        else if(aa == "Trp")
            pos_is_good = PositionIsTrp(seq, i);
        else
            VERIFY_MSG(false, "Unknown aa identifier: " << aa);
        if(pos_is_good)
            aa_positions.push_back(i);
    }
    return aa_positions;
}

bool SimpleCdr3Calculator::Cdr3PositionsCanBeRefined(IgIsotype isotype,
                                                     std::pair<size_t, size_t> new_pos,
                                                     std::pair<size_t, size_t> old_pos) {
    size_t avg_hc_cdr3_length = 45;
    size_t avg_lc_cdr3_length = 27;
    size_t avg_cdr3_length = avg_hc_cdr3_length;
    if(isotype.IsLightChain())
        avg_cdr3_length = avg_lc_cdr3_length;
    size_t old_len = old_pos.second - old_pos.first + 1;
    size_t new_len = new_pos.second - new_pos.first + 1;
    return (new_pos.first > old_pos.first) and
        (abs_diff(new_len, avg_cdr3_length) <= abs_diff(old_len, avg_cdr3_length));
}

std::string SimpleCdr3Calculator::FindCdr3Positions(const IsotypeUmiSequence& umi_record) {
    auto cys_pos = ComputeAaPositions(umi_record.sequence, "Cys");
    auto phe_pos = ComputeAaPositions(umi_record.sequence, "Phe");
    auto trp_pos = ComputeAaPositions(umi_record.sequence, "Trp");
    std::pair<size_t, size_t> cdr3_pos(0, 0);
    for(auto cit = cys_pos.begin(); cit != cys_pos.end(); cit++) {
        size_t start_pos = *cit + 3;
        for(auto pit = phe_pos.begin(); pit != phe_pos.end(); pit++) {
            size_t end_pos = *pit - 1;
            if(start_pos >= end_pos)
                continue;
            if(Cdr3PositionsCanBeRefined(umi_record.isotype,
                                         std::make_pair(start_pos, end_pos),
                                         cdr3_pos))
                cdr3_pos = std::make_pair(start_pos, end_pos);
        }
        for(auto tit = trp_pos.begin(); tit != trp_pos.end(); tit++) {
            size_t end_pos = *tit - 1;
            if(start_pos >= end_pos)
                continue;
            if(Cdr3PositionsCanBeRefined(umi_record.isotype,
                                         std::make_pair(start_pos, end_pos),
                                         cdr3_pos))
                cdr3_pos = std::make_pair(start_pos, end_pos);
        }
    }

    //INFO("CDR3 positions for " << umi_record.sequence << " are (" <<
    //             cdr3_pos.first << ", " <<
    //             cdr3_pos.second << ")");
    //auto cdr3 = seqan::infixWithLength(umi_record.sequence,
    //                                   cdr3_pos.first,
    //                                   cdr3_pos.second - cdr3_pos.first + 1);
    //INFO("CDR3 sequence: " << cdr3);
    seqan::String<char> cdr3_char;
    seqan::assign(cdr3_char, umi_record.sequence);
    return std::string(seqan::toCString(cdr3_char)).substr(cdr3_pos.first,
                                                                     cdr3_pos.second - cdr3_pos.first + 1);
}
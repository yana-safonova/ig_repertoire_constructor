#include <verify.hpp>

#include <boost/algorithm/string.hpp>

#include "amino_acid_motif.hpp"

namespace core {
    std::string AminoAcidExtendedAlphabet = "ACDEFGHIKLMNPQRSTVWYX";

    bool CharIsExtendedAminoAcid(char sym) {
        for(auto it = AminoAcidExtendedAlphabet.begin(); it != AminoAcidExtendedAlphabet.end(); it++)
            if(*it == sym)
                return true;
        return false;
    }

    bool StringIsExtendedAminoAcid(std::string aa_seq) {
        for(auto it = aa_seq.begin(); it != aa_seq.end(); it++)
            if(!CharIsExtendedAminoAcid(*it))
                return false;
        return true;
    }

    void AminoAcidMotifs::Initialize() {
        std::vector<std::string> alternative_splits;
        boost::split(alternative_splits, aa_motif_, boost::is_any_of("/"));
        for(auto it = alternative_splits.begin(); it != alternative_splits.end(); it++) {
            std::cout << *it << std::endl;
            std::vector<std::string> seq_splits;
            boost::split(seq_splits, *it, boost::is_any_of("-"));
            std::string aa_motif;
            for(auto it2 = seq_splits.begin(); it2 != seq_splits.end(); it2++)
                aa_motif += *it2;
            VERIFY_MSG(StringIsExtendedAminoAcid(aa_motif), "String " << aa_motif << " is not extended aa sequence");
            aa_variants_.push_back(aa_motif);
        }
    }
}
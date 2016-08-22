#pragma once

#include <vector>
#include <string>

namespace core {
    // X means unknown amino acid
    extern std::string AminoAcidExtendedAlphabet;

    bool CharIsExtendedAminoAcid(char sym);

    bool StringIsExtendedAminoAcid(std::string aa_seq);

    class AminoAcidMotifs {
        std::string aa_motif_;
        std::vector<std::string> aa_variants_;

        void Initialize();

    public:
        AminoAcidMotifs(std::string aa_motif) : aa_motif_(aa_motif) {
            Initialize();
        }

        size_t size() const { return aa_variants_.size(); }

        typedef std::vector<std::string>::const_iterator AaVariantsConstIterator;

        AaVariantsConstIterator cbegin() const { return aa_variants_.begin(); }

        AaVariantsConstIterator cend() const { return aa_variants_.end(); }

        // todo: refactor it!
        size_t length() const { return aa_variants_[0].size(); }
    };
}
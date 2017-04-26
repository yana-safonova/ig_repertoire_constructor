//
// Created by Andrew Bzikadze on 4/25/17.
//

#pragma once

#include <array>
#include <algorithm>
#include <string>
#include <annotation_utils/cdr_labeling_primitives.hpp>

namespace ig_simulator {

class FastStopCodonCheckerDetails {
friend class FastStopCodonChecker;
private:
    constexpr static unsigned get_hash(const char *s,
                                       const unsigned hash_base,
                                       const unsigned hash_base_sq) {
        return s[0] + hash_base * s[1] + hash_base_sq * s[2];
    }

    static unsigned get_hash(std::string &&s,
                             const unsigned hash_base,
                             const unsigned hash_base_sq) {
        return get_hash(s.c_str(), hash_base, hash_base_sq);
    }
};

class FastStopCodonChecker {
private:
    constexpr const static unsigned hash_base { 10 };
    constexpr const static unsigned hash_base_sq { hash_base * hash_base };

    /**
     *  I have to use the following hack and refrain from using std::array for stop_codons because
     *  on OSX El Capitan neigher std::array::operator[] nor std::array std::get
     *  are not declared constexpr.
     *  The following code should be used instead in the future:
     *  @code
     *  constexpr const static std::array<const char[4], 3> stop_codons { "TAG", "TAA", "TGA" };
     *  constexpr const static std::array<unsigned, 3> stop_codons_hashes
     *      {{
     *           sc_checker_details::get_hash(std::get<0>(stop_codons), hash_base, hash_base_sq),
     *           sc_checker_details::get_hash(std::get<1>(stop_codons), hash_base, hash_base_sq),
     *           sc_checker_details::get_hash(std::get<2>(stop_codons), hash_base, hash_base_sq)
     *      }};
     */
    constexpr const static std::array<unsigned, 3> stop_codons_hashes
        {{
             FastStopCodonCheckerDetails::get_hash("TAG", hash_base, hash_base_sq),
             FastStopCodonCheckerDetails::get_hash("TAA", hash_base, hash_base_sq),
             FastStopCodonCheckerDetails::get_hash("TGA", hash_base, hash_base_sq)
        }};

public:
    bool static HasStopCodon(const std::string& str, size_t orf) {
        for(size_t i = orf; i + 2 < str.length(); i += 3) {
            size_t hash = FastStopCodonCheckerDetails::get_hash(str.substr(i, 3),
                                                       hash_base, hash_base_sq);
            if (std::find(stop_codons_hashes.begin(), stop_codons_hashes.end(), hash)
                != stop_codons_hashes.end())
            {
                return true;
            }
        }
        return false;
    }

    bool static HasStopCodon(const std::string& str, const annotation_utils::CDRLabeling& labeling) {
        return HasStopCodon(str, labeling.cdr1.start_pos % 3);
    }

    FastStopCodonChecker() = delete;
    FastStopCodonChecker(const FastStopCodonChecker&) = delete;
    FastStopCodonChecker(FastStopCodonChecker&&) = delete;
    FastStopCodonChecker& operator=(const FastStopCodonChecker&) = delete;
    FastStopCodonChecker& operator=(FastStopCodonChecker&&) = delete;
};

} // End namespace ig_simulator

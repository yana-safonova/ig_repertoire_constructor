#pragma once

#include <bitset>
#include <cassert>
#include <seqan/seq_io.h>
#include <unordered_map>
#include <vector>

namespace fast_ig_tools {

template <size_t MAX_LEN = 640>
class bitread {
public:
    template <typename TString>
    explicit bitread(const TString &s) : len__{seqan::length(s)} {
        assert(len__ <= MAX_LEN);

        for (size_t i = 0; i < len__; i++) {
            std::bitset<2> tmp = seqan::ordValue(s[i]);
            shifted0__.evenbit[i] = tmp[0];
            shifted0__.oddbit[i] = tmp[1];
        }
    }

    bitread() : bitread("") {}

    size_t length() const {
        return len__;
    }

    // Thread-safe
    size_t dist(const bitread &s) const {
        size_t mlen = std::min<size_t>(s.length(), this->length());
        const auto &shifted = shifted0__;
        std::bitset<MAX_LEN> res = ((s.shifted0__.evenbit ^ shifted.evenbit) | (s.shifted0__.oddbit ^ shifted.oddbit)) & this->right_masks__(MAX_LEN - mlen);
        return res.count();
    }

    // Thread-safe
    size_t dist_equal_len(const bitread &s) const {
        assert(this->length() == s.length());
        const auto &shifted = shifted0__;
        std::bitset<MAX_LEN> res = ((s.shifted0__.evenbit ^ shifted.evenbit) | (s.shifted0__.oddbit ^ shifted.oddbit));
        return res.count();
    }

    struct BitReadPair__ {
        std::bitset<MAX_LEN> oddbit;
        std::bitset<MAX_LEN> evenbit;

        BitReadPair__ operator<<(size_t i) const {
            BitReadPair__ tmp;
            tmp.evenbit = evenbit << i;
            tmp.oddbit = oddbit << i;
            return tmp;
        }

        BitReadPair__ operator>>(size_t i) const {
            BitReadPair__ tmp;
            tmp.evenbit = evenbit >> i;
            tmp.oddbit = oddbit >> i;
            return tmp;
        }

        BitReadPair__ shifted(int i) const {
            if (i >= 0) {
                return (*this) >> i;
            } else {
                return (*this) << (-i);
            }
        }
    };

    const BitReadPair__& content() const {
        return shifted0__;
    }

protected:
    static const std::bitset<MAX_LEN> *get_masks_right__() {
        static std::bitset<MAX_LEN> right_masks[MAX_LEN + 1];

        right_masks[0].set();
        for (size_t i = 1; i < MAX_LEN; i++) {
            right_masks[i] = right_masks[i - 1] >> 1;
        }

        return right_masks;
    }

    static const std::bitset<MAX_LEN> *get_masks_left__() {
        static std::bitset<MAX_LEN> left_masks[MAX_LEN + 1];

        left_masks[0].set();
        for (size_t i = 1; i < MAX_LEN; i++) {
            left_masks[i] = left_masks[i - 1] << 1;
        }

        return left_masks;
    }

    static const std::bitset<MAX_LEN> &right_masks__(size_t i) {
        static auto m = get_masks_right__();
        return m[i];
    }

    static const std::bitset<MAX_LEN> &left_masks__(size_t i) {
        static auto m = get_masks_left__();
        return m[i];
    }

private:
    size_t len__;
    BitReadPair__ shifted0__;
};

template <size_t MAX_LEN = 640>
class shifted_bitread : public bitread<MAX_LEN> {
public:
    using Base = bitread<MAX_LEN>;
    // Inherit cinstructor
    using Base::bitread;

    // Not thread-safe!!!
    size_t dist(const Base &s, int shift) {
        assert(-shift <= static_cast<int>(this->length()));

        size_t mlen = std::min<size_t>(s.length(), this->length() - shift);
        const auto &shifted = getShifted__(shift);

        std::bitset<MAX_LEN> res = ((s.content().evenbit ^ shifted.evenbit) | (s.content().oddbit ^ shifted.oddbit)) & this->right_masks__(MAX_LEN - mlen);
        if (shift < 0) {
            res &= this->left_masks__(-shift);
        }
        return res.count();
    }

private:
    using BitReadPair__ = typename Base::BitReadPair__;

    const BitReadPair__ &getShifted__(int shift) {
        auto it = shifted_reads__.find(shift);

        if (it == shifted_reads__.cend()) {
            return shifted_reads__[shift] = Base::content().shifted(shift);
        } else {
            return it->second;
        }
    }

    std::unordered_map<int, BitReadPair__> shifted_reads__;
};

}  // namespace fast_ig_tools

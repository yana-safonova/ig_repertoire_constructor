#pragma once

#include <vector>
#include <verify.hpp>

namespace algorithms {

    struct KmerMatch {
        int needle_pos;
        int read_pos;

        int shift() const {
            return read_pos - needle_pos;
        }

        static bool less_shift(const KmerMatch &m1, const KmerMatch &m2) {
            return (m1.shift() < m2.shift()) || (m1.shift() == m2.shift() && m1.read_pos < m2.read_pos);
        }
    };

    class SubjectKmerMatches {
        std::vector<std::vector<KmerMatch>> kmer_matches_;
        size_t num_subjects_;

    public:
        SubjectKmerMatches(size_t num_subjects) :
                kmer_matches_(num_subjects),
                num_subjects_(num_subjects) { }

        void Update(size_t subject_index, KmerMatch kmer_match) {
            VERIFY(subject_index < num_subjects_);
            kmer_matches_[subject_index].push_back(kmer_match);
        }

        size_t size() const { return kmer_matches_.size(); }

        std::vector<KmerMatch>& operator[](size_t subject_index) {
            VERIFY(subject_index < num_subjects_);
            return kmer_matches_[subject_index];
        }

        typedef std::vector<std::vector<KmerMatch>>::const_iterator SubjectKmerMatchesIterator;

        SubjectKmerMatchesIterator cbegin() const { return kmer_matches_.cbegin(); }

        SubjectKmerMatchesIterator cend() const { return kmer_matches_.cend(); }
    };

    template<typename SubjectDatabase, typename StringType>
    class KmerIndexHelper {
    protected:
        const SubjectDatabase &db_;

    public:
        KmerIndexHelper(const SubjectDatabase &db) : db_(db) { }

        virtual StringType GetDbRecordByIndex(size_t index) const = 0;// {
        //    return db_[index];
        //}

        virtual size_t GetStringLength(const StringType &s) const = 0; // {
        //    return seqan::length(s);
        //}

        virtual size_t GetDbSize() const = 0;
    };
}
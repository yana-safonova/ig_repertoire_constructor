#pragma once

#include <read_archive.hpp>
#include <germline_utils/germline_databases/immune_gene_database.hpp>
#include <block_alignment/pairwise_block_alignment.hpp>
#include <germline_utils/germline_databases/custom_gene_database.hpp>

#include <boost/optional.hpp>

namespace vj_finder {
    class ImmuneGeneHit {
    protected:
        const core::Read *read_ptr_;
        const germline_utils::ImmuneGene* immune_gene_ptr_;
        algorithms::PairwiseBlockAlignment block_alignment_;
        bool strand_;

    public:
        ImmuneGeneHit() :
            read_ptr_(NULL),
            immune_gene_ptr_(NULL),
            block_alignment_(),
            strand_(false) { }

        ImmuneGeneHit(const core::Read& read,
                      const germline_utils::ImmuneGene &immune_gene,
                      algorithms::PairwiseBlockAlignment &block_alignment,
                      bool strand) :
                read_ptr_(&read),
                immune_gene_ptr_(&immune_gene),
                block_alignment_(std::move(block_alignment)),
                strand_(strand) { }

        bool Strand() const { return strand_; }

        const core::Read& Read() const {
            VERIFY(!Empty());
            return *read_ptr_;
        }

        const algorithms::PairwiseBlockAlignment& BlockAlignment() const {
            return block_alignment_;
        }

        const germline_utils::ImmuneGene& ImmuneGene() const {
            VERIFY(!Empty());
            return *immune_gene_ptr_;
        }

        germline_utils::ChainType Chain()  const {
            VERIFY(!Empty());
            return immune_gene_ptr_->Chain();
        }

        virtual size_t SegmentLength() const = 0;

        virtual int LeftUncovered() const { return 0; }

        virtual int RightUncovered() const { return 0; }

        double Score() const {
            VERIFY(!Empty());
            return block_alignment_.score;
        }

        virtual int Start() const = 0;

        virtual int End() const = 0;

        bool Empty() const { return read_ptr_ == NULL or immune_gene_ptr_ == NULL; }


        // first match positions are inclusive
        size_t FirstMatchReadPos() const { return block_alignment_.first_match_read_pos(); }

        size_t FirstMatchGenePos() const { return block_alignment_.first_match_subject_pos(); }

        // last match positions are exclusive
        size_t LastMatchReadPos() const { return block_alignment_.last_match_read_pos(); }

        size_t LastMatchGenePos() const { return block_alignment_.last_match_subject_pos(); }


        virtual void AddShift(int shift) {
            block_alignment_.add_read_shift(shift);
        }

        virtual void ExtendFirstMatch(int left_shift) {
            block_alignment_.extend_first_match(left_shift);
        }

        virtual void ExtendLastMatch(int right_shift) {
            block_alignment_.extend_last_match(right_shift);
        }
    };

    class VGeneHit : public ImmuneGeneHit {
    public:
        VGeneHit(const core::Read& read,
                 const germline_utils::ImmuneGene &immune_gene,
                 algorithms::PairwiseBlockAlignment &block_alignment,
                 bool strand) : ImmuneGeneHit(read, immune_gene, block_alignment, strand) { }

        virtual size_t SegmentLength() const {
            VERIFY(!Empty());
            return block_alignment_.left_half_segment_length();
        }

        virtual int LeftUncovered() const {
            VERIFY(!Empty());
            return std::max(0, -block_alignment_.start());
        }

        virtual int Start() const {
            VERIFY(!Empty());
            return block_alignment_.start();
        }

        virtual int End() const {
            VERIFY(!Empty());
            return int(block_alignment_.last_match_read_pos());
        }
    };

    class JGeneHit : public ImmuneGeneHit {
    public:
        JGeneHit(const core::Read& read,
                 const germline_utils::ImmuneGene &immune_gene,
                 algorithms::PairwiseBlockAlignment &block_alignment,
                 bool strand) : ImmuneGeneHit(read, immune_gene, block_alignment, strand) { }

        virtual size_t SegmentLength() const {
            VERIFY(!Empty());
            return block_alignment_.right_half_segment_length();
        }

        virtual int RightUncovered() const {
            VERIFY(!Empty());
            //std::cout << block_alignment_.finish() << " - " << read_.length() << std::endl;
            return std::max(0, block_alignment_.finish() - static_cast<int>(read_ptr_->length()));
        }

        virtual int Start() const {
            VERIFY(!Empty());
            return int(block_alignment_.first_match_read_pos());
        }

        virtual int End() const {
            VERIFY(!Empty());
            return block_alignment_.finish();
        }
    };

    class VJHits {
        const core::Read *read_ptr_;
        std::vector<VGeneHit> v_hits_;
        std::vector<JGeneHit> j_hits_;

        void CheckConsistencyFatal(const VGeneHit &v_hit);

        void CheckConsistencyFatal(const JGeneHit &j_hit);

    public:
        VJHits(const core::Read &read) : read_ptr_(&read) { }

        void AddVHit(VGeneHit v_hit) {
            v_hits_.push_back(v_hit);
        }

        void AddJHit(JGeneHit j_hit) {
            j_hits_.push_back(j_hit);
        }

        size_t NumVHits() const { return v_hits_.size(); }

        size_t NumJHits() const { return j_hits_.size(); }

        VGeneHit GetVHitByIndex(size_t index) const {
            VERIFY_MSG(index < NumVHits(), "Index " << index << " exceeds number V hits");
            return v_hits_[index];
        }

        JGeneHit GetJHitByIndex(size_t index) const {
            VERIFY_MSG(index < NumJHits(), "Index " << index << " exceeds number J hits");
            return j_hits_[index];
        }

        size_t AlignedSegmentLength() const {
            return GetJHitByIndex(0).End() - GetVHitByIndex(0).Start();
        }

        const core::Read& Read() const { return *read_ptr_; }

        void AddLeftShift(int shift) {
            for(auto it = v_hits_.begin(); it != v_hits_.end(); it++)
                it->AddShift(shift);
        }

        void AddRightShift(int shift) {
            for(auto it = j_hits_.begin(); it != j_hits_.end(); it++)
                it->AddShift(shift);
        }

        void ExtendFirstMatch(int left_shift) {
            for(auto it = v_hits_.begin(); it != v_hits_.end(); it++) {
                it->ExtendFirstMatch(left_shift);
            }
        }

        void ExtendLastMatch(int right_shift) {
            for(auto it = j_hits_.begin(); it != j_hits_.end(); it++) {
                it->ExtendLastMatch(right_shift);
            }
        }
    };

    class CustomGermlineDbHelper : public algorithms::KmerIndexHelper<germline_utils::CustomGeneDatabase,
            seqan::Dna5String> {
    public:
        CustomGermlineDbHelper(const germline_utils::CustomGeneDatabase& db) :
                KmerIndexHelper<germline_utils::CustomGeneDatabase, seqan::Dna5String>(db) { }

        seqan::Dna5String GetDbRecordByIndex(size_t index) const {
            return db_[index].seq();
        }

        size_t GetStringLength(const seqan::Dna5String &s) const {
            return seqan::length(s);
        }

        size_t GetDbSize() const {
            return db_.size();
        }
    };

    class ImmuneGeneGermlineDbHelper : public algorithms::KmerIndexHelper<germline_utils::ImmuneGeneDatabase,
            seqan::Dna5String> {
    public:
        ImmuneGeneGermlineDbHelper(const germline_utils::ImmuneGeneDatabase& db) :
                KmerIndexHelper<germline_utils::ImmuneGeneDatabase, seqan::Dna5String>(db) { }

        seqan::Dna5String GetDbRecordByIndex(size_t index) const {
            return db_[index].seq();
        }

        size_t GetStringLength(const seqan::Dna5String &s) const {
            return seqan::length(s);
        }

        size_t GetDbSize() const {
            return db_.size();
        }
    };
}
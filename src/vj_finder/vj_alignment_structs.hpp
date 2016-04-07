#pragma once

#include <read_archive.hpp>
#include <germline_utils/germline_databases/immune_gene_database.hpp>
#include <block_alignment/pairwise_block_alignment.hpp>
#include <germline_utils/germline_databases/custom_gene_database.hpp>

namespace vj_finder {
    class ImmuneGeneHit {
    protected:
        const core::Read* read_ptr_;
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

        virtual size_t LeftUncovered() const { return 0; }

        virtual size_t RightUncovered() const { return 0; }

        double Score() const {
            VERIFY(!Empty());
            return block_alignment_.score;
        }

        virtual int Start() const = 0;

        virtual int End() const = 0;

        bool Empty() const { return read_ptr_ == NULL or immune_gene_ptr_ == NULL; }
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

        virtual size_t LeftUncovered() const {
            VERIFY(!Empty());
            return std::max(0, block_alignment_.start);
        }

        virtual int Start() const {
            VERIFY(!Empty());
            return block_alignment_.start;
        }

        virtual int End() const {
            VERIFY(!Empty());
            return block_alignment_.last_match_read_pos();
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

        virtual size_t RightUncovered() const {
            VERIFY(!Empty());
            return std::max(0, block_alignment_.finish - static_cast<int>(read_ptr_->length()));
        }

        virtual int Start() const {
            VERIFY(!Empty());
            return block_alignment_.first_match_read_pos();
        }

        virtual int End() const {
            VERIFY(!Empty());
            return block_alignment_.finish;
        }
    };

    class VJHit {
        const core::Read& read_;
        VGeneHit v_hit_;
        JGeneHit j_hit_;

        void CheckConsistency();

    public:
        VJHit(const core::Read& read,
              VGeneHit v_hit,
              JGeneHit j_hit) :
                read_(read),
                v_hit_(v_hit),
                j_hit_(j_hit) { }

        const core::Read& Read() const { return read_; }

        bool Strand() const { return v_hit_.Strand(); }

        germline_utils::ChainType Chain() const { return v_hit_.Chain(); }

        const ImmuneGeneHit& V() const { return v_hit_; }

        const ImmuneGeneHit& J() const { return j_hit_; }

        int FinalLength() const {
            return j_hit_.End() - v_hit_.Start();
        }
    };

    class VJHits {
        const core::Read &read_;
        std::vector<VGeneHit> v_hits_;
        std::vector<JGeneHit> j_hits_;

        void CheckConsistencyFatal(const VGeneHit &v_hit);

        void CheckConsistencyFatal(const JGeneHit &j_hit);

    public:
        VJHits(const core::Read &read) : read_(read) { }

        void AddVHit(VGeneHit v_hit) {
            v_hits_.push_back(v_hit);
        }

        void AddJHit(JGeneHit j_hit) {
            j_hits_.push_back(j_hit);
        }

        size_t NumVHits() { return v_hits_.size(); }

        size_t NumJHits() { return j_hits_.size(); }

        //const VJHit& operator[](size_t index) {
        //    VERIFY_MSG(index < size(), "Index " << index << " exceeds VJ hits size");
        //    return vj_hits_[index];
        //}
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

/*
    class GermlineDbReadKmerHashes : public algorithms::SubjectQueryKmerIndex<germline_utils::CustomGeneDatabase,
            seqan::Dna5String, core::Read> {
    public:
        GermlineDbReadKmerHashes(const germline_utils::CustomGeneDatabase &db,
                                 size_t k) :
                algorithms::SubjectQueryKmerHashes<germline_utils::CustomGeneDatabase,
                        seqan::Dna5String, core::Read>(db, k) { }

    protected:
        virtual void Initialize() {
            size_t gene_index = 0;
            for(auto it = db_.cbegin(); it != db_.cend(); it++) {
                auto immune_gene_db = db_.GetDbByGeneType(*it);
                for(auto it2 = immune_gene_db.cbegin(); it2 != immune_gene_db.cend(); it2++) {
                    auto gene_string = it2->seq();
                    auto hashes = algorithms::polyhashes(gene_string, k_);
                    UpdateMap(hashes, GetStringLength(gene_string), gene_index);
                    gene_index++;
                }
            }
        }

        virtual seqan::Dna5String GetQueryString(const core::Read &q) const {
            return q.seq;
        }
    };
    */

    /*
    class VJAlignmentOld {
        Dna5String read;
        char strand;
        std::string locus;
        std::vector <BlockAligner::Alignment> v_hits, j_hits;
        std::shared_ptr<const GermlineLociVJDB> vbase, jbase;

    public:
        VJAlignment(const Dna5String &read,
                    char strand,
                    const std::string &locus,
                    const std::vector <BlockAligner::Alignment> &v_hits,
                    const std::vector <BlockAligner::Alignment> &j_hits,
                    std::shared_ptr<const GermlineLociVJDB> vbase,
                    std::shared_ptr<const GermlineLociVJDB> jbase) : read{read},
                                                                     strand{strand},
                                                                     locus{locus},
                                                                     v_hits{v_hits},
                                                                     j_hits{j_hits},
                                                                     vbase{vbase},
                                                                     jbase{jbase} {
            ;;;
        }

        static const VJAlignment &EmptyVJAlignment() {
            static VJAlignment empty;

            return empty;
        }

        const Dna5String &Read() const {
            return read;
        }

        char Strand() const {
            return strand;
        }

        const std::string &Locus() const {
            return locus;
        }

        bool empty() const {
            return v_hits.empty() || j_hits.empty();
        }

        operator bool() const {
            return !empty();
        }

        size_t VHitsSize() const {
            return v_hits.size();
        }

        size_t JHitsSize() const {
            return j_hits.size();
        }

        const BlockAligner::Alignment &VHit(size_t v_hit_index = 0) const {
            assert(v_hit_index < v_hits.size());

            return v_hits[v_hit_index];
        }

        const BlockAligner::Alignment &JHit(size_t j_hit_index = 0) const {
            assert(j_hit_index < j_hits.size());

            return j_hits[j_hit_index];
        }

        template<typename Tparam>
        Dna5String FixCropFill(const Tparam &param,
                               size_t v_hit_index = 0, size_t j_hit_index = 0) const {
            return FixCropFill(param.fix_left, param.crop_left, param.fill_left,
                               param.fix_right, param.crop_right, param.fill_right,
                               v_hit_index, j_hit_index);
        }

        Dna5String FixCropFill(size_t fix_left, bool crop_left, bool fill_left,
                               size_t fix_right, bool crop_right, bool fill_right,
                               size_t v_hit_index = 0, size_t j_hit_index = 0) const {
            const auto &v_hit = VHit(v_hit_index);
            const auto &j_hit = JHit(j_hit_index);

            const auto &v_gene = VSeq(v_hit_index);
            const auto &j_gene = JSeq(j_hit_index);

            int left_shift = v_hit.path.left_shift();
            int right_shift = j_hit.path.right_shift();

            Dna5String result = read;

            for (size_t i = 0; i < fix_left; ++i) {
                int gene_pos = static_cast<int>(i) - left_shift;

                if (gene_pos >= 0 && gene_pos < static_cast<int>(length(v_gene))) {
                    result[i] = v_gene[gene_pos];
                }
            }

            for (size_t i = length(read) - fix_right; i < length(read); ++i) {
                int gene_pos = static_cast<int>(i) - right_shift;

                if (gene_pos >= 0 && gene_pos < static_cast<int>(length(j_gene))) {
                    result[i] = j_gene[gene_pos];
                }
            }

            if (crop_right && j_hit.finish < static_cast<int>(length(result))) {
                result = seqan::prefix(result, j_hit.finish);
            } else if (fill_right && j_hit.finish > static_cast<int>(length(result))) {
                result += seqan::suffix(j_gene, length(j_gene) - (j_hit.finish - length(result)));
            }

            if (crop_left && v_hit.start > 0) {
                result = seqan::suffix(result, v_hit.start);
            } else if (fill_left && v_hit.start < 0) {
                Dna5String pre = seqan::prefix(v_gene, -v_hit.start);
                pre += result;
                result = pre;
            }

            return result;
        }

        const Dna5String &VSeq(size_t v_hit_index = 0) const {
            const auto &v_hit = VHit(v_hit_index);
            return vbase->v_reads[v_hit.needle_index];
        }

        const Dna5String &JSeq(size_t j_hit_index = 0) const {
            const auto &j_hit = JHit(j_hit_index);
            return jbase->j_reads[j_hit.needle_index];
        }

        size_t VSegmentLength(size_t v_hit_index = 0) const {
            return VHit(v_hit_index).left_half_segment_length();
        }

        size_t JSegmentLength(size_t j_hit_index = 0) const {
            return JHit(j_hit_index).right_half_segment_length();
        }

        size_t LeftUncovered(size_t v_hit_index = 0) const {
            return std::max(0, -VHit(v_hit_index).start);
        }

        size_t RightUncovered(size_t j_hit_index = 0) const {
            return std::max(0, JHit(j_hit_index).finish - static_cast<int>(length(read)));
        }

        double VScore(size_t v_hit_index = 0) const {
            return VHit(v_hit_index).score;
        }

        double JScore(size_t j_hit_index = 0) const {
            return JHit(j_hit_index).score;
        }

        int VStart(size_t v_hit_index = 0) const {
            return VHit(v_hit_index).start;
        }

        int JStart(size_t j_hit_index = 0) const {
            return JHit(j_hit_index).first_match_read_pos();
        }

        int VEnd(size_t v_hit_index = 0) const {
            return VHit(v_hit_index).last_match_read_pos();
        }

        int JEnd(size_t j_hit_index = 0) const {
            return JHit(j_hit_index).finish;
        }

        size_t FinalLength(size_t v_hit_index = 0,
                           size_t j_hit_index = 0) const {
            return JEnd(j_hit_index) - VStart(v_hit_index);
        }

        TAlign VAlignmentSeqAn(size_t v_hit_index = 0) const {
            return VHit(v_hit_index).seqan_alignment(read, VSeq(v_hit_index));
        }

        TAlign JAlignmentSeqAn(size_t j_hit_index = 0) const {
            return JHit(j_hit_index).seqan_alignment(read, JSeq(j_hit_index));
        }

        std::string VMatches(size_t v_hit_index = 0) const {
            return VHit(v_hit_index).visualize(read, VSeq(v_hit_index));
        }

        std::string JMatches(size_t j_hit_index = 0) const {
            return JHit(j_hit_index).visualize(read, JSeq(j_hit_index));
        }


        std::pair <TAlign, TAlign> AlignmentsSeqan(size_t v_hit_index = 0,
                                                   size_t j_hit_index = 0) const {
            return {VAlignmentSeqAn(v_hit_index), JAlignmentSeqAn(j_hit_index)};
        }

        const CharString &VId(size_t v_hit_index = 0) const {
            return vbase->v_ids[VHit(v_hit_index).needle_index];
        }

        const CharString &JId(size_t j_hit_index = 0) const {
            return jbase->j_ids[JHit(j_hit_index).needle_index];
        }

    private:
        VJAlignment() = default;
    };
     */
}
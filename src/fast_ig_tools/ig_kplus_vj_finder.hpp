#pragma once

#include <cassert>
#include <stdexcept>
#include <vector>
#include <memory>

using std::vector;

#include "fast_ig_tools.hpp"

#include <seqan/seq_io.h>
using seqan::Dna5String;
using seqan::CharString;
using seqan::length;

#include "ig_block_alignment.hpp"

namespace fast_ig_tools {

struct GermlineLociVJDB {
    std::vector<Dna5String> v_reads;
    std::vector<CharString> v_ids;
    std::vector<Dna5String> j_reads;
    std::vector<CharString> j_ids;

    GermlineLociVJDB& extend(const GermlineLociVJDB &db) {
        v_reads.insert(v_reads.end(), db.v_reads.cbegin(), db.v_reads.cend());
        v_ids.insert(v_ids.end(), db.v_ids.cbegin(), db.v_ids.cend());

        j_reads.insert(j_reads.end(), db.j_reads.cbegin(), db.j_reads.cend());
        j_ids.insert(j_ids.end(), db.j_ids.cbegin(), db.j_ids.cend());

        return *this;
    }
};


struct VJAlignment {
    size_t locus_id;
    std::string locus;

    char strand;
    Dna5String read;
    const GermlineLociVJDB *vbase, *jbase;

    std::vector<BlockAligner::Alignment> v_hits, j_hits;

    const Dna5String& Read() const {
        return read;
    }

    char Strand() const {
        return strand;
    }

    const std::string& Locus() const {
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

    const BlockAligner::Alignment& VHit(size_t v_hit_index = 0) const {
        assert(v_hit_index < v_hits.size());

        return v_hits[v_hit_index];
    }

    const BlockAligner::Alignment& JHit(size_t j_hit_index = 0) const {
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

    const Dna5String& VSeq(size_t v_hit_index = 0) const {
        const auto &v_hit = VHit(v_hit_index);
        return vbase->v_reads[v_hit.needle_index];
    }

    const Dna5String& JSeq(size_t j_hit_index = 0) const {
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

    std::pair<TAlign, TAlign> AlignmentsSeqan(size_t v_hit_index = 0,
                                              size_t j_hit_index = 0) const {
        return { VAlignmentSeqAn(v_hit_index), JAlignmentSeqAn(j_hit_index) };
    }

    const CharString& VId(size_t v_hit_index = 0) const {
        return vbase->v_ids[VHit(v_hit_index).needle_index];
    }

    const CharString& JId(size_t j_hit_index = 0) const {
        return jbase->j_ids[JHit(j_hit_index).needle_index];
    }
};


class VJAligner {
public:
    template<typename Tparam>
    VJAligner(const Tparam &param) : db_directory{param.db_directory},
                                     organism{param.organism},
                                     pseudogenes{param.pseudogenes},
                                     locus_names{expand_loci({ param.loci })} {
        for (const std::string &locus : locus_names) {
            GermlineLociVJDB db;

            {
                std::string v_file = gene_file_name(locus, "V", pseudogenes);
                seqan::SeqFileIn seqFileIn(v_file.c_str());
                readRecords(db.v_ids, db.v_reads, seqFileIn);
            }

            {
                std::string j_file = gene_file_name(locus, "J", pseudogenes);
                seqan::SeqFileIn seqFileIn(j_file.c_str());
                readRecords(db.j_ids, db.j_reads, seqFileIn);
            }

            locus_databases.push_back(db);
        }

        for (size_t i = 0; i < locus_databases.size(); ++i) {
            const auto &db = locus_databases[i];

            all_loci_database.extend(db);

            for (size_t _ = 0; _ < db.v_reads.size(); ++_) {
                locus_index.push_back(i);
            }
        }

        valigner.reset(new BlockAligner(all_loci_database.v_reads, param.K, param.max_global_gap,
                                        100500, 100500,
                                        param.max_local_insertions, param.max_local_deletions, param.min_k_coverage));
        jaligner.reset(new BlockAligner(all_loci_database.j_reads, param.word_size_j, param.max_global_gap,
                                        100500, 100500,
                                        param.max_local_insertions, param.max_local_deletions, param.min_k_coverage_j));

        for (const auto db : locus_databases) {
            BlockAligner *p = new BlockAligner(db.j_reads, param.word_size_j, param.max_global_gap,
                                               100500, 100500,
                                               param.max_local_insertions, param.max_local_deletions, param.min_k_coverage_j);
            jaligners.push_back(std::shared_ptr<BlockAligner>(p));
        }
    }

    size_t vbase_size() const {
        return all_loci_database.v_reads.size();
    }

    size_t jbase_size() const {
        return all_loci_database.j_reads.size();
    }


    template<typename Tparam>
    VJAlignment Query(const Dna5String &read,
                      const Tparam &param) const {
        bool fix_strand = param.fix_strand;
        bool consistent_loci = param.consistent_loci;

        VJAlignment result_alignment;

        std::vector<BlockAligner::Alignment> result;
        Dna5String stranded_read;
        char strand;

        // Fix strand if asked
        std::tie(result, stranded_read, strand) = fix_strand ? correct_strand(read, param.max_candidates) : correct_strand_fake(read, param.max_candidates);

        if (result.empty()) {
            return result_alignment; // Empty TODO Add marker
        }


        result_alignment.read = stranded_read;
        result_alignment.strand = strand;
        result_alignment.vbase = &all_loci_database;
        // Aling V --- already done


        // Identify locus, in case of misconsensus report empty alignment
        std::vector<size_t> locus_ids;
        locus_ids.reserve(result.size());
        for (const auto &align : result) {
            size_t loc_ind = locus_index[align.needle_index];

            locus_ids.push_back(loc_ind);
        }

        // Check equality
        if (consistent_loci && !all_equal(locus_ids)) {
            return result_alignment; // Empty TODO Add marker
        }

        // Save V
        result_alignment.v_hits = result; // FIXME use move()

        size_t locus_id = result_alignment.locus_id = locus_ids[0];
        result_alignment.locus = locus_names[locus_id];
        result_alignment.jbase = &locus_databases[locus_id];

        // Aling V --- already done
        // Find minimum suffix after V gene
        int end_of_v = max_map(result.cbegin(), result.cend(),
                               [](const BlockAligner::Alignment &align) -> int { return align.last_match_read_pos(); });

        // Align J
        // auto jaligns = jaligners[locus_id]->query(seqan::suffix(read, end_of_v), param.max_candidates_j);
        const auto &jal = consistent_loci ? jaligners[locus_id] : jaligner;
        auto jaligns = jal->query(stranded_read, param.max_candidates_j, end_of_v);

        result_alignment.j_hits = jaligns;

        // Report
        return result_alignment;
    }

private:
    template<typename T>
    static bool all_equal(const T &v) {
        for (size_t i = 1; i < v.size(); ++i) {
            if (v[i] != v[0]) {
                return false;
            }
        }

        return true;
    }


    template<typename Titer, typename Tf>
    static auto max_map(Titer b, Titer e, const Tf &f) -> decltype(f(*b)) { // TODO Add decay
        assert(b != e);

        auto result = f(*b++);
        for(;b != e; ++b) {
            auto cur = f(*b);
            if (result <  cur) {
                result = cur;
            }
        }

        return result;
    }


    std::tuple<std::vector<BlockAligner::Alignment>, Dna5String, char> correct_strand(const Dna5String &read,
                                                                                      size_t max_candidates) const {
        Dna5String read_rc = read;
        reverseComplement(read_rc);

        auto result_pstrand = valigner->query(read, max_candidates);
        auto result_nstrand = valigner->query(read_rc, max_candidates);

        int pscore = (result_pstrand.size() > 0) ? result_pstrand[0].kp_coverage : 0;
        int nscore = (result_nstrand.size() > 0) ? result_nstrand[0].kp_coverage : 0;

        char strand  = (pscore >= nscore) ? '+' : '-';
        const auto &stranded_read = (strand == '+') ? read : read_rc;
        const auto &result = (strand == '+') ? result_pstrand : result_nstrand;

        return std::make_tuple(result, stranded_read, strand);
    }


    std::tuple<std::vector<BlockAligner::Alignment>, Dna5String, char> correct_strand_fake(const Dna5String &read,
                                                                                           size_t max_candidates) const {
        auto result_pstrand = valigner->query(read, max_candidates);

        return std::make_tuple(result_pstrand, read, '+' );
    }


    static std::vector<std::string> expand_loci(const std::vector<std::string> &loci) {
        std::vector<std::string> IG = { "IGH", "IGK", "IGL" };
        std::vector<std::string> TR = { "TRA", "TRB", "TRD", "TRG" };
        std::vector<std::string> ALL(IG);
        ALL.insert(ALL.end(), TR.cbegin(), TR.cend());

        std::vector<std::string> result;
        for (const std::string &locus : loci) {
            if (locus == "IG") {
                result.insert(result.end(), IG.cbegin(), IG.cend());
            } else if (locus == "TR") {
                result.insert(result.end(), TR.cbegin(), TR.cend());
            } else if (locus == "all") {
                result.insert(result.end(), ALL.cbegin(), ALL.cend());
            } else if (std::find(ALL.cbegin(), ALL.cend(), locus) != ALL.cend()) {
                result.push_back(locus);
            } else {
                throw std::invalid_argument("Invalid locus name: " + locus);
            }
        }

        // Remove duplicates
        remove_duplicates(result);

        return result;
    }


    std::string gene_file_name(const std::string &locus,
                               const std::string &gene,
                               bool pseudo = false) const {
        return db_directory + "/" + organism + "/" + locus + gene + (pseudo ? "-allP" : "") + ".fa";
    }


    std::vector<GermlineLociVJDB> locus_databases;
    GermlineLociVJDB all_loci_database;

    const std::string db_directory;
    const std::string organism;
    bool pseudogenes;
    std::vector<std::string> locus_names;

    std::vector<size_t> locus_index;

    std::shared_ptr<BlockAligner> valigner;
    std::shared_ptr<BlockAligner> jaligner;
    std::vector<std::shared_ptr<BlockAligner>> jaligners;
};

} // namespace fast_ig_tools

// vim: ts=4:sw=4

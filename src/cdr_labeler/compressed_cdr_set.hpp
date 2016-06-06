#pragma once

#include <vj_alignment_info.hpp>
#include <convert.hpp>
#include <annotation_utils/annotated_clone_set.hpp>

namespace cdr_labeler {
    struct CDRKey {
        seqan::CharString v_name;
        seqan::CharString j_name;
        seqan::Dna5String cdr_seq;
        size_t id;

        CDRKey(seqan::CharString v_name, seqan::CharString j_name,
               seqan::Dna5String cdr_seq, size_t id) : v_name(v_name),
                                                       j_name(j_name),
                                                       cdr_seq(cdr_seq),
                                                       id(id) { }

        bool operator==(const CDRKey& obj) const {
            return v_name == obj.v_name and j_name == obj.j_name and cdr_seq == obj.cdr_seq;
        }
    };

    struct CDRKeyHasher {
        size_t operator()(const CDRKey &obj) const {
            return std::hash<std::string>()(std::string(seqan::toCString(obj.v_name))) *
                    std::hash<std::string>()(std::string(seqan::toCString(obj.j_name))) *
                    std::hash<std::string>()(core::dna5String_to_string(obj.cdr_seq));
        }
    };

    class CompressedCDRSet {
        annotation_utils::StructuralRegion region_;
        const vj_finder::VJAlignmentInfo &alignment_info_;
        const annotation_utils::CDRAnnotatedCloneSet &clone_set_;

        std::unordered_map<CDRKey, size_t, CDRKeyHasher> compressed_cdrs_map_;

    public:
        typedef std::pair<CDRKey, size_t> CompressedCDR;

    private:
        std::vector<CompressedCDR> compressed_cdrs_;

        void Initialize();

    public:
        CompressedCDRSet(annotation_utils::StructuralRegion region,
                         const vj_finder::VJAlignmentInfo &alignment_info,
                         const annotation_utils::CDRAnnotatedCloneSet &clone_set) : region_(region),
                                                                                    alignment_info_(alignment_info),
                                                                                    clone_set_(clone_set) {
            Initialize();
        }

        typedef std::vector<CompressedCDR>::const_iterator CompressedCDRConstIterator;

        CompressedCDRConstIterator cbegin() const { return compressed_cdrs_.cbegin(); }

        CompressedCDRConstIterator cend() const { return compressed_cdrs_.cend(); }

        size_t size() const { return compressed_cdrs_.size(); }
    };
}
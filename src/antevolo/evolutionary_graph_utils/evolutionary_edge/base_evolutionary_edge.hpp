#pragma once

#include <annotation_utils/annotated_clone.hpp>
namespace  antevolo {

    template<typename T>
    size_t HammingDistance(T seq1, T seq2) {
        size_t dist = 0;
        size_t min_length = std::min<size_t>(seqan::length(seq1), seqan::length(seq2));
        for(size_t i = 0; i < min_length; i++) {
            if(seq1[i] != seq2[i])
                dist++;
        }
        return dist;
    }

    enum EvolutionaryEdgeType {
        UnknownEdgeType,
        DirectedEdgeType,
        UndirectedEdgeType,
        IntersectedEdgeType,
        ReverseDirectedEdgeType
    };

    class BaseEvolutionaryEdge {
    protected:
        EvolutionaryEdgeType edge_type;
        std::shared_ptr<annotation_utils::AnnotatedClone> src_clone;
        std::shared_ptr<annotation_utils::AnnotatedClone> dst_clone;
        size_t src_clone_num;
        size_t dst_clone_num;
        bool synonymous;
        size_t cdr3_distance;

    public:
        BaseEvolutionaryEdge(const annotation_utils::AnnotatedClone &src_clone,
                             const annotation_utils::AnnotatedClone &dst_clone,
                             size_t src_num, size_t dst_num)
                : edge_type(EvolutionaryEdgeType::UnknownEdgeType),
                  src_clone(std::make_shared<annotation_utils::AnnotatedClone>(src_clone)),
                  dst_clone(std::make_shared<annotation_utils::AnnotatedClone>(dst_clone)),
                  src_clone_num(src_num),
                  dst_clone_num(dst_num),
                  cdr3_distance(0) {

            synonymous = true;
            cdr3_distance = HammingDistance(this->src_clone->CDR3(), this->dst_clone->CDR3());
            auto src_AA_seq = this->src_clone->AA();
            auto dst_AA_seq = this->dst_clone->AA();
            size_t src_length = seqan::length(src_AA_seq);
            if (src_length != seqan::length(dst_AA_seq)) {
                synonymous = false;
            }
            for (size_t i = 0; i < src_length; ++i) {
                if (src_AA_seq[i] != dst_AA_seq[i]) {
                    synonymous = false;
                }
            }

        }

        virtual const std::shared_ptr<annotation_utils::AnnotatedClone>& SrcClone() const {
            return src_clone;
        };
        virtual const std::shared_ptr<annotation_utils::AnnotatedClone>& DstClone() const {
            return dst_clone;
        };
        virtual size_t SrcNum() const {
            return src_clone_num;
        }
        virtual size_t DstNum() const {
            return dst_clone_num;
        }

        virtual bool IsSynonymous() const { return synonymous; }

        virtual size_t CDR3Distance() const { return cdr3_distance; }

        virtual size_t NumAddedShms() const { return size_t(-1); }

        virtual size_t NumSharedShms() const { return size_t(-1); }

        virtual bool Empty() const { return edge_type == EvolutionaryEdgeType::UnknownEdgeType; }

        virtual bool IsDirected() const { return false; }

        virtual bool IsUndirected() const { return false; }

        virtual bool IsIntersected() const { return false; }

        virtual bool IsDoubleMutated() const { return false; }

        virtual bool IsReverseDirected() const { return false; }

        virtual std::string TypeString() const { return "unknown"; }

        virtual size_t Length() const { return size_t(-1); }

        virtual ~BaseEvolutionaryEdge() {}

//        bool operator==(const EvolutionaryEdge &edge) const {
//            return src_clone == edge.src_clone and dst_clone == edge.dst_clone and edge_type == edge.edge_type;
//        }
    };

    typedef std::shared_ptr<BaseEvolutionaryEdge> EvolutionaryEdgePtr;
}


#include <vj_query_aligner.hpp>
#include "annotated_clone_by_read_constructor.hpp"

namespace antevolo {
    annotation_utils::AnnotatedClone AnnotatedCloneByReadConstructor::GetCloneByRead(core::Read& read) const {

        vj_finder::VJQueryAligner vj_query_aligner(vj_finder_params_, labeled_v_db_, labeled_j_db_);
        vj_finder::VJHits vj_hits = vj_query_aligner.Align(read);
        cdr_labeler::ReadCDRLabeler read_labeler(shm_config_, v_labeling_, j_labeling_);
        return read_labeler.CreateAnnotatedClone(vj_hits);
    }
    annotation_utils::AnnotatedClone AnnotatedCloneByReadConstructor::GetCloneByReadWithSpecificGenes (
            core::Read &read,
            const germline_utils::ImmuneGene &v_gene,
            const germline_utils::ImmuneGene &j_gene) const {
        germline_utils::CustomGeneDatabase single_V_gene_db(germline_utils::SegmentType::VariableSegment);
        single_V_gene_db.AddImmuneGene(v_gene);
        germline_utils::CustomGeneDatabase single_J_gene_db(germline_utils::SegmentType::JoinSegment);
        single_J_gene_db.AddImmuneGene(j_gene);

        vj_finder::VJQueryAligner vj_query_aligner(vj_finder_params_, single_V_gene_db, single_J_gene_db);
        vj_finder::VJHits pre_vj_hits = vj_query_aligner.Align(read);

        vj_finder::VJHits vj_hits(pre_vj_hits.Read());
        for (size_t i = 0; i < pre_vj_hits.NumVHits(); ++i) {
            auto pre_v_hit = pre_vj_hits.GetVHitByIndex(i);
            const auto& gene_type = pre_v_hit.ImmuneGene().GeneType();
            auto& immune_gene_db = labeled_v_db_.GetDbByGeneType(gene_type);
            size_t gene_index = immune_gene_db.GetIndexByName(pre_v_hit.ImmuneGene().name());
            algorithms::PairwiseBlockAlignment block_alignment = pre_v_hit.BlockAlignment();
            auto v_hit = vj_finder::VGeneHit(vj_hits.Read(),
                                             immune_gene_db.GetImmuneGeneByIndex(gene_index),
                                             block_alignment,
                                             pre_v_hit.Strand());
            vj_hits.AddVHit(v_hit);
        }
        for (size_t i = 0; i < pre_vj_hits.NumJHits(); ++i) {
            auto pre_j_hit = pre_vj_hits.GetJHitByIndex(i);
            const auto& gene_type = pre_j_hit.ImmuneGene().GeneType();
            auto& immune_gene_db = labeled_j_db_.GetDbByGeneType(gene_type);
            size_t gene_index = immune_gene_db.GetIndexByName(pre_j_hit.ImmuneGene().name());
            algorithms::PairwiseBlockAlignment block_alignment = pre_j_hit.BlockAlignment();
            auto j_hit = vj_finder::JGeneHit(vj_hits.Read(),
                                             immune_gene_db.GetImmuneGeneByIndex(gene_index),
                                             block_alignment,
                                             pre_j_hit.Strand());
            vj_hits.AddJHit(j_hit);
        }
        cdr_labeler::ReadCDRLabeler read_labeler(shm_config_, v_labeling_, j_labeling_);
        return read_labeler.CreateAnnotatedClone(vj_hits);
    }
    annotation_utils::AnnotatedClone AnnotatedCloneByReadConstructor::GetCloneByReadAndAlignment(
            std::tuple<core::Read,
                    seqan::Align<seqan::Dna5String>,
                    seqan::Align<seqan::Dna5String>> tpl,
            const germline_utils::ImmuneGene &v_gene,
            const germline_utils::ImmuneGene &j_gene) const {
        auto& read = std::get<0>(tpl);
        auto& seqan_v_alignment = std::get<1>(tpl);
        auto& seqan_j_alignment = std::get<2>(tpl);
        alignment_utils::ImmuneGeneReadAlignment v_alignment(&v_gene,
                                                             &read,
                                                             seqan_v_alignment,
                                                             0);
        alignment_utils::ImmuneGeneReadAlignment j_alignment(&j_gene,
                                                             &read,
                                                             seqan_j_alignment,
                                                             0);
        auto v_cdr_labeling = v_labeling_.GetLabelingByGene(v_gene);
        annotation_utils::CDRRange read_cdr1(v_alignment.QueryPositionBySubjectPosition(v_cdr_labeling.cdr1.start_pos),
                                             v_alignment.QueryPositionBySubjectPosition(v_cdr_labeling.cdr1.end_pos));
        annotation_utils::CDRRange read_cdr2(v_alignment.QueryPositionBySubjectPosition(v_cdr_labeling.cdr2.start_pos),
                                             v_alignment.QueryPositionBySubjectPosition(v_cdr_labeling.cdr2.end_pos));
        auto j_cdr_labeling = j_labeling_.GetLabelingByGene(j_gene);
        annotation_utils::CDRRange read_cdr3(v_alignment.QueryPositionBySubjectPosition(v_cdr_labeling.cdr3.start_pos),
                                             j_alignment.QueryPositionBySubjectPosition(j_cdr_labeling.cdr3.end_pos));
        annotation_utils::CDRLabeling cdr_labeling(read_cdr1, read_cdr2, read_cdr3);
        cdr_labeler::ReadCDRLabeler read_labeler(shm_config_, v_labeling_, j_labeling_);
        return read_labeler.GetCloneCalculator().ComputeAnnotatedClone(read,
                                                                 cdr_labeling,
                                                                 v_alignment,
                                                                 j_alignment);

    }

    std::pair<size_t, size_t> AnnotatedCloneByReadConstructor::GetGeneCDR3Ranges(
            const germline_utils::ImmuneGene& v_gene,
            const germline_utils::ImmuneGene& j_gene) const {
        return {v_labeling_.GetLabelingByGene(v_gene).cdr3.start_pos,
                j_labeling_.GetLabelingByGene(j_gene).cdr3.end_pos};
    }
}

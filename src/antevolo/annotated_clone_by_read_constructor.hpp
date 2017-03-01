#pragma once


#include <annotation_utils/annotated_clone.hpp>
#include <read_labeler.hpp>
#include <vj_finder_config.hpp>

namespace antevolo {

    class AnnotatedCloneByReadConstructor {
        germline_utils::CustomGeneDatabase& labeled_v_db_;
        germline_utils::CustomGeneDatabase& labeled_j_db_;
        const cdr_labeler::DbCDRLabeling& v_labeling_;
        const cdr_labeler::DbCDRLabeling& j_labeling_;
        const vj_finder::VJFinderConfig::AlgorithmParams& vj_finder_params_;
        const cdr_labeler::CDRLabelerConfig::SHMFindingParams &shm_config_;
    public:
        AnnotatedCloneByReadConstructor(germline_utils::CustomGeneDatabase& labeled_v_db,
                                        germline_utils::CustomGeneDatabase& labeled_j_db,
                                        const cdr_labeler::DbCDRLabeling& v_labeling,
                                        const cdr_labeler::DbCDRLabeling& j_labeling,
                                        const vj_finder::VJFinderConfig::AlgorithmParams& vj_finder_params,
                                        const cdr_labeler::CDRLabelerConfig::SHMFindingParams &shm_config) :
                labeled_v_db_(labeled_v_db),
                labeled_j_db_(labeled_j_db),
                v_labeling_(v_labeling),
                j_labeling_(j_labeling),
                vj_finder_params_(vj_finder_params),
                shm_config_(shm_config) {}

        annotation_utils::AnnotatedClone GetCloneByRead(core::Read& read) const;

        annotation_utils::AnnotatedClone GetCloneByReadWithSpecificGenes(
                core::Read &read,
                const germline_utils::ImmuneGene &v_gene,
                const germline_utils::ImmuneGene &j_gene) const;

        annotation_utils::AnnotatedClone GetCloneByReadAndAlignment(
                std::tuple<core::Read,
                        seqan::Align<seqan::Dna5String>,
                        seqan::Align<seqan::Dna5String>> tpl,
                const germline_utils::ImmuneGene &v_gene,
                const germline_utils::ImmuneGene &j_gene) const;

        const vj_finder::VJFinderConfig::AlgorithmParams& GetVJFinderParams() const {
            return vj_finder_params_;
        }

        std::pair<size_t, size_t> GetGeneCDR3Ranges(const germline_utils::ImmuneGene& v_gene,
                                                    const germline_utils::ImmuneGene& j_gene) const;

    };

}

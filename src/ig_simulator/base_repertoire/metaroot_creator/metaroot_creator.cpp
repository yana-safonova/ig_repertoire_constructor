//
// Created by Andrew Bzikadze on 3/20/17.
//

#include "metaroot_creator.hpp"

#include <random>

#include <base_repertoire/gene_chooser/abstract_gene_chooser.hpp>
#include "random_generator.hpp"
#include "annotation_utils/cdr_labeling_primitives.hpp"

#include <iostream>

namespace ig_simulator {

AbstractMetarootCreator::AbstractMetarootCreator(const MetarootSimulationParams& config,
                                                 std::vector<germline_utils::CustomGeneDatabase>& db,
                                                 AbstractVDJGeneChooserCPtr&& gene_chooser):
    v_db_p(&db.front()),
    j_db_p(&db.back()),
    prob_cleavage_v(check_probability(config.cleavage_params.prob_cleavage_v)),
    prob_cleavage_j(check_probability(config.cleavage_params.prob_cleavage_j)),
    gene_chooser_p(std::move(gene_chooser)),
    nucl_remover_p(get_nucleotides_remover(config.nucleotides_remover_params)),
    nucl_creator_p(get_nucleotides_creator(config.p_nucleotides_creator_params)),
    nucl_inserter_p(get_nucleotides_inserter(config.n_nucleotides_inserter_params)),
    v_cdr_db(cdr_labeler::GermlineDbLabeler(db.front(), config.cdr_labeler_config.cdrs_params).ComputeLabeling()),
    j_cdr_db(cdr_labeler::GermlineDbLabeler(db.back(),  config.cdr_labeler_config.cdrs_params).ComputeLabeling()),
    productivity_checker()
{
    VERIFY(db.size() >= 2);
    VERIFY(v_db_p->size() > 0);
    VERIFY(j_db_p->size() > 0);
}

VJMetarootCreator::VJMetarootCreator(const MetarootSimulationParams& config,
                                     std::vector<germline_utils::CustomGeneDatabase>& db):
    AbstractMetarootCreator(config, db, get_gene_chooser(config.gene_chooser_params, db))
{
    VERIFY(db.size() == 2);
}

AbstractMetarootCPtr VJMetarootCreator::Createroot() const {
    auto genes_ind = gene_chooser_p->ChooseGenes();
    VERIFY(std::get<1>(genes_ind) == size_t(-1));

    bool is_cleavage_v       = std::bernoulli_distribution(prob_cleavage_v)(MTSingleton::GetInstance());
    bool is_cleavage_j       = std::bernoulli_distribution(prob_cleavage_j)(MTSingleton::GetInstance());

    int cleavage_v = is_cleavage_v ?
                     static_cast<int>(nucl_remover_p->RemoveInVGene()) :
                     -static_cast<int>(nucl_creator_p->CreateInVGene());

    int cleavage_j = is_cleavage_j ?
                     static_cast<int>(nucl_remover_p->RemoveInJGene()) :
                     -static_cast<int>(nucl_creator_p->CreateInJGene());

    seqan::Dna5String vj_insertion(nucl_inserter_p->GetVJInsertion());

    const auto& v_gene = (*v_db_p)[std::get<0>(genes_ind)];
    const auto& j_gene = (*j_db_p)[std::get<2>(genes_ind)];

    annotation_utils::CDRLabeling cdr_labeling(v_cdr_db.GetLabelingByGene(v_gene));
    annotation_utils::CDRLabeling j_gene_cdr_labeling(j_cdr_db.GetLabelingByGene(j_gene));

    if (not v_cdr_db.CDRLabelingIsEmpty(v_gene) and not j_cdr_db.CDRLabelingIsEmpty(j_gene)) {
        long long cdr3_end = static_cast<long long>(v_gene.length()) +
            -cleavage_v +
            static_cast<long long>(seqan::length(vj_insertion)) +
            -cleavage_j +
            static_cast<long long>(j_gene_cdr_labeling.cdr3.end_pos);
        VERIFY(cdr3_end >= 0);
        cdr_labeling.cdr3.end_pos = static_cast<size_t>(cdr3_end);
    }

    VJMetaroot metaroot { v_db_p, j_db_p,
                          std::get<0>(genes_ind), std::get<2>(genes_ind),
                          cdr_labeling,
                          cleavage_v, cleavage_j,
                          vj_insertion };
    if (not productivity_checker.IsProductive(metaroot)) {
        metaroot.SetNonProductive();
    }
    return AbstractMetarootCPtr(new VJMetaroot(std::move(metaroot)));
}

VDJMetarootCreator::VDJMetarootCreator(const MetarootSimulationParams& config,
                                       std::vector<germline_utils::CustomGeneDatabase>& db):
    AbstractMetarootCreator(config, db, get_gene_chooser(config.gene_chooser_params, db)),
    d_db_p(&db.at(1)),
    prob_cleavage_d_left(check_probability(config.cleavage_params.prob_cleavage_d_left)),
    prob_cleavage_d_right(check_probability(config.cleavage_params.prob_cleavage_d_right))
{
    VERIFY(db.size() == 3);
    VERIFY(d_db_p->size() > 0);
}

AbstractMetarootCPtr VDJMetarootCreator::Createroot() const {
    auto genes_ind = gene_chooser_p->ChooseGenes();

    bool is_cleavage_v       = std::bernoulli_distribution(prob_cleavage_v)(MTSingleton::GetInstance());
    bool is_cleavage_d_left  = std::bernoulli_distribution(prob_cleavage_d_left)(MTSingleton::GetInstance());
    bool is_cleavage_d_right = std::bernoulli_distribution(prob_cleavage_d_right)(MTSingleton::GetInstance());
    bool is_cleavage_j       = std::bernoulli_distribution(prob_cleavage_j)(MTSingleton::GetInstance());

    int cleavage_v = is_cleavage_v ?
                     static_cast<int>(nucl_remover_p->RemoveInVGene()) :
                     -static_cast<int>(nucl_creator_p->CreateInVGene());

    int cleavage_d_left = is_cleavage_d_left ?
                          static_cast<int>(nucl_remover_p->RemoveInDGeneLeft()) :
                          -static_cast<int>(nucl_creator_p->CreateInDGeneLeft());

    int cleavage_d_right = is_cleavage_d_right ?
                          static_cast<int>(nucl_remover_p->RemoveInDGeneRight()) :
                          -static_cast<int>(nucl_creator_p->CreateInDGeneRight());

    int cleavage_j = is_cleavage_j ?
                     static_cast<int>(nucl_remover_p->RemoveInJGene()) :
                     -static_cast<int>(nucl_creator_p->CreateInJGene());

    seqan::Dna5String vd_insertion(nucl_inserter_p->GetVDInsertion());
    seqan::Dna5String dj_insertion(nucl_inserter_p->GetDJInsertion());

    const auto& v_gene = (*v_db_p)[std::get<0>(genes_ind)];
    const auto& d_gene = (*d_db_p)[std::get<1>(genes_ind)];
    const auto& j_gene = (*j_db_p)[std::get<2>(genes_ind)];

    annotation_utils::CDRLabeling cdr_labeling(v_cdr_db.GetLabelingByGene(v_gene));
    annotation_utils::CDRLabeling j_gene_cdr_labeling(j_cdr_db.GetLabelingByGene(j_gene));

    if (not v_cdr_db.CDRLabelingIsEmpty(v_gene) and not j_cdr_db.CDRLabelingIsEmpty(j_gene)) {
        long long cdr3_end = static_cast<long long>(v_gene.length()) +
            -cleavage_v +
            static_cast<long long>(seqan::length(vd_insertion)) +
            -cleavage_d_left +
            static_cast<long long>(d_gene.length()) +
            -cleavage_d_right +
            static_cast<long long>(seqan::length(dj_insertion)) +
            -cleavage_j +
            static_cast<long long>(j_gene_cdr_labeling.cdr3.end_pos);

        VERIFY(cdr3_end >= 0);
        cdr_labeling.cdr3.end_pos = static_cast<size_t>(cdr3_end);
    }

    VDJMetaroot metaroot { v_db_p, d_db_p, j_db_p,
                           std::get<0>(genes_ind), std::get<1>(genes_ind), std::get<2>(genes_ind),
                           cdr_labeling,
                           cleavage_v, cleavage_d_left, cleavage_d_right, cleavage_j,
                           vd_insertion, dj_insertion };

    if (not productivity_checker.IsProductive(metaroot)) {
        metaroot.SetNonProductive();
    }
    return AbstractMetarootCPtr(new VDJMetaroot(std::move(metaroot)));
}

} // End namespace ig_simulator

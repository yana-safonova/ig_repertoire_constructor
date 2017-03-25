//
// Created by Andrew Bzikadze on 3/20/17.
//

#include "metaroot_creator.hpp"

#include <random>

#include <gene_chooser/abstract_gene_chooser.hpp>
#include "random_generator.hpp"

#include <iostream>

namespace ig_simulator {

VJMetaRoot VJMetarootCreator::CreateRoot() const {
    auto genes_ind = gene_chooser_p->ChooseGenes();
    VERIFY(std::get<1>(genes_ind) == size_t(-1));

    bool is_cleavage_v       = std::bernoulli_distribution(prob_cleavage_v)(MTSingleton::GetInstance());
    bool is_cleavage_j       = std::bernoulli_distribution(prob_cleavage_j)(MTSingleton::GetInstance());

    int cleavage_v = is_cleavage_v ?
                     -static_cast<int>(nucl_remover_p->RemoveInVGene()) :
                     static_cast<int>(nucl_creator_p->CreateInVGene());

    int cleavage_j = is_cleavage_j ?
                     -static_cast<int>(nucl_remover_p->RemoveInJGene()) :
                     static_cast<int>(nucl_creator_p->CreateInJGene());

    return VJMetaRoot(v_db_p, j_db_p,
                      std::get<0>(genes_ind), std::get<2>(genes_ind),
                      cleavage_v, cleavage_j);
}

VDJMetaRoot VDJMetarootCreator::CreateRoot() const {
    auto genes_ind = gene_chooser_p->ChooseGenes();

    bool is_cleavage_v       = std::bernoulli_distribution(prob_cleavage_v)(MTSingleton::GetInstance());
    bool is_cleavage_d_left  = std::bernoulli_distribution(prob_cleavage_d_left)(MTSingleton::GetInstance());
    bool is_cleavage_d_right = std::bernoulli_distribution(prob_cleavage_d_right)(MTSingleton::GetInstance());
    bool is_cleavage_j       = std::bernoulli_distribution(prob_cleavage_j)(MTSingleton::GetInstance());

    int cleavage_v = is_cleavage_v ?
                     -static_cast<int>(nucl_remover_p->RemoveInVGene()) :
                     static_cast<int>(nucl_creator_p->CreateInVGene());

    int cleavage_d_left = is_cleavage_d_left ?
                          -static_cast<int>(nucl_remover_p->RemoveInDGeneLeft()) :
                          static_cast<int>(nucl_creator_p->CreateInDGeneLeft());

    int cleavage_d_right = is_cleavage_d_right ?
                          -static_cast<int>(nucl_remover_p->RemoveInDGeneRight()) :
                          static_cast<int>(nucl_creator_p->CreateInDGeneRight());

    int cleavage_j = is_cleavage_j ?
                     -static_cast<int>(nucl_remover_p->RemoveInJGene()) :
                     static_cast<int>(nucl_creator_p->CreateInJGene());

    return VDJMetaRoot(v_db_p, d_db_p, j_db_p,
                       std::get<0>(genes_ind), std::get<1>(genes_ind), std::get<2>(genes_ind),
                       cleavage_v, cleavage_d_left, cleavage_d_right, cleavage_j,
                       nucl_inserter_p->GetVDInsertion(), nucl_inserter_p->GetDJInsertion());

}

} // End namespace ig_simulator

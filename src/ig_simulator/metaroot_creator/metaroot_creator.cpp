//
// Created by Andrew Bzikadze on 3/20/17.
//

#include "metaroot_creator.hpp"

#include <random>

#include <gene_chooser/abstract_gene_chooser.hpp>
#include "random_generator.hpp"

#include <iostream>

namespace ig_simulator {

void MetarootCreator::CreateRoot() const {
    auto genes_ind = gene_chooser_p->ChooseGenes();
    std::cout << std::get<0>(genes_ind) << " " << std::get<1>(genes_ind) << " " << std::get<2>(genes_ind) << std::endl;

    std::bernoulli_distribution bern_d(prob_cleavage);
    bool is_cleavage = bern_d(MTSingleton::GetInstance());
    std::cout << is_cleavage << std::endl;

    std::cout << nucl_remover_p->RemoveInVGene() << std::endl;
}

} // End namespace ig_simulator

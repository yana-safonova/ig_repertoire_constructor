//
// Created by Andrew Bzikadze on 3/20/17.
//

#pragma once

#include "germline_utils/chain_type.hpp"
#include "base_repertoire/gene_chooser/abstract_gene_chooser.hpp"
#include "base_repertoire/gene_chooser/config_based_getter.hpp"
#include "base_repertoire/nucleotides_remover/abstract_nucleotides_remover.hpp"
#include "base_repertoire/nucleotides_remover/config_based_getter.hpp"
#include "base_repertoire/p_nucleotides_creator/abstract_nucleotides_creator.hpp"
#include "base_repertoire/p_nucleotides_creator/config_based_getter.hpp"
#include "base_repertoire/n_nucleotides_inserter/abstract_n_nucleotides_inserter.hpp"
#include "base_repertoire/n_nucleotides_inserter/config_based_getter.hpp"
#include "base_repertoire/metaroot/metaroot.hpp"
#include "germline_db_labeler.hpp"
#include "germline_db_labeling.hpp"
#include "cdr_config.hpp"
#include "base_repertoire/productivity_checker/productivity_checker.hpp"

namespace ig_simulator {

class AbstractMetarootCreator {
protected:
    // TODO
    // Databases are not declared `const` since cdr_labeler requires: see germline_db_labeler.hpp
    // @code: DbCDRLabeling GermlineDbLabeler::ComputeLabeling();
    // is not declared const
    germline_utils::CustomGeneDatabase * v_db_p;
    germline_utils::CustomGeneDatabase * j_db_p;

    const double prob_cleavage_v;
    const double prob_cleavage_j;

    const AbstractVDJGeneChooserCPtr gene_chooser_p;
    const AbstractNucleotidesRemoverCPtr nucl_remover_p;
    const AbstractPNucleotidesCreatorCPtr nucl_creator_p;
    const AbstractNNucleotidesInserterCPtr nucl_inserter_p;

    const cdr_labeler::DbCDRLabeling v_cdr_db;
    const cdr_labeler::DbCDRLabeling j_cdr_db;

    const ProductivityChecker productivity_checker;

public:
    AbstractMetarootCreator(const MetarootSimulationParams& config,
                            std::vector<germline_utils::CustomGeneDatabase>& db,
                            AbstractVDJGeneChooserCPtr&& gene_chooser);

    AbstractMetarootCreator() = delete;
    AbstractMetarootCreator(const AbstractMetarootCreator&) = delete;
    AbstractMetarootCreator(AbstractMetarootCreator&&) = delete;
    AbstractMetarootCreator& operator=(const AbstractMetarootCreator&) = delete;
    AbstractMetarootCreator& operator=(AbstractMetarootCreator&) = delete;

    virtual AbstractMetarootCPtr Createroot() const = 0;
    virtual ~AbstractMetarootCreator() { }
};
using AbstractMetarootCreatorCPtr = std::unique_ptr<const AbstractMetarootCreator>;


class VJMetarootCreator final : public AbstractMetarootCreator {
public:

    VJMetarootCreator(const MetarootSimulationParams& config,
                      std::vector<germline_utils::CustomGeneDatabase>& db);

    AbstractMetarootCPtr Createroot() const override;
};


class VDJMetarootCreator final : public AbstractMetarootCreator {
private:
    germline_utils::CustomGeneDatabase * d_db_p;

    const double prob_cleavage_d_left;
    const double prob_cleavage_d_right;

public:
    VDJMetarootCreator(const MetarootSimulationParams& config,
                       std::vector<germline_utils::CustomGeneDatabase>& db);

    AbstractMetarootCPtr Createroot() const override;
};

} // End namespace ig_simulator

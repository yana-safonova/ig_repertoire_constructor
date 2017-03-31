//
// Created by Andrew Bzikadze on 3/20/17.
//

#pragma once

#include "germline_utils/chain_type.hpp"
#include "gene_chooser/abstract_gene_chooser.hpp"
#include "gene_chooser/config_based_getter.hpp"
#include "nucleotides_remover/abstract_nucleotides_remover.hpp"
#include "nucleotides_remover/config_based_getter.hpp"
#include "p_nucleotides_creator/abstract_nucleotides_creator.hpp"
#include "p_nucleotides_creator/config_based_getter.hpp"
#include "n_nucleotides_inserter/abstract_n_nucleotides_inserter.hpp"
#include "n_nucleotides_inserter/config_based_getter.hpp"
#include "metaroot/metaroot.hpp"
#include "germline_db_labeler.hpp"
#include "germline_db_labeling.hpp"
#include "cdr_config.hpp"

namespace ig_simulator {

class AbstractMetarootCreator {
protected:
    // TODO
    // Databases are not declared `const` since cdr_labeler requires: see germline_db_labeler.hpp
    // @code: DbCDRLabeling GermlineDbLabeler::ComputeLabeling();
    // is not declared const
    germline_utils::CustomGeneDatabase * v_db_p;
    germline_utils::CustomGeneDatabase * j_db_p;

    double prob_cleavage_v;
    double prob_cleavage_j;

    AbstractVDJGeneChooserPtr gene_chooser_p;
    AbstractNucleotidesRemoverPtr nucl_remover_p;
    AbstractPNucleotidesCreatorPtr nucl_creator_p;
    AbstractNNucleotidesInserterPtr nucl_inserter_p;

    cdr_labeler::DbCDRLabeling v_cdr_db;
    cdr_labeler::DbCDRLabeling j_cdr_db;

public:
    AbstractMetarootCreator(const MetarootSimulationParams& config,
                            germline_utils::CustomGeneDatabase *v_db_p,
                            germline_utils::CustomGeneDatabase *j_db_p,
                            AbstractVDJGeneChooserPtr&& gene_chooser):
        v_db_p(check_pointer(v_db_p)),
        j_db_p(check_pointer(j_db_p)),
        prob_cleavage_v(config.cleavage_params.prob_cleavage_v),
        prob_cleavage_j(config.cleavage_params.prob_cleavage_j),
        gene_chooser_p(std::move(gene_chooser)),
        nucl_remover_p(get_nucleotides_remover(config.nucleotides_remover_params)),
        nucl_creator_p(get_nucleotides_creator(config.p_nucleotides_creator_params)),
        nucl_inserter_p(get_nucleotides_inserter(config.n_nucleotides_inserter_params)),
        v_cdr_db(cdr_labeler::GermlineDbLabeler(*check_pointer(v_db_p), config.cdr_labeler_config.cdrs_params).ComputeLabeling()),
        j_cdr_db(cdr_labeler::GermlineDbLabeler(*check_pointer(j_db_p), config.cdr_labeler_config.cdrs_params).ComputeLabeling())
    {
        VERIFY(v_db_p->size() > 0);
        VERIFY(j_db_p->size() > 0);
        VERIFY(prob_cleavage_v >= 0 and prob_cleavage_v <= 1);
        VERIFY(prob_cleavage_j >= 0 and prob_cleavage_j <= 1);
    }

    virtual AbstractMetarootPtr Createroot() const = 0;
    virtual ~AbstractMetarootCreator() { }
};

class VJMetarootCreator : public AbstractMetarootCreator {
public:

    VJMetarootCreator(const MetarootSimulationParams& config,
                      germline_utils::CustomGeneDatabase *v_db_p,
                      germline_utils::CustomGeneDatabase *j_db_p):
        AbstractMetarootCreator(config, v_db_p, j_db_p, get_gene_chooser(config.gene_chooser_params, {v_db_p, j_db_p}))
    { }

    virtual AbstractMetarootPtr Createroot() const override;
};

class VDJMetarootCreator : public AbstractMetarootCreator {
private:
    germline_utils::CustomGeneDatabase * d_db_p;

    const double prob_cleavage_d_left;
    const double prob_cleavage_d_right;

public:
    VDJMetarootCreator(const MetarootSimulationParams& config,
                       germline_utils::CustomGeneDatabase *v_db_p,
                       germline_utils::CustomGeneDatabase *j_db_p,
                       germline_utils::CustomGeneDatabase *d_db_p):
        AbstractMetarootCreator(config, v_db_p, j_db_p, get_gene_chooser(config, {v_db_p, d_db_p, j_db_p})),
        d_db_p(check_pointer(d_db_p)),
        prob_cleavage_d_left(config.cleavage_params.prob_cleavage_d_left),
        prob_cleavage_d_right(config.cleavage_params.prob_cleavage_d_right)
    {
        VERIFY(d_db_p->size() > 0);
        VERIFY(prob_cleavage_d_left >= 0 and prob_cleavage_d_left <= 1);
        VERIFY(prob_cleavage_d_right >= 0 and prob_cleavage_d_right <= 1);
    }

    virtual AbstractMetarootPtr Createroot() const override;
};

} // End namespace ig_simulator

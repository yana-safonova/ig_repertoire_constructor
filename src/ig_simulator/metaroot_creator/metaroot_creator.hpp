//
// Created by Andrew Bzikadze on 3/20/17.
//

#pragma once

#include "germline_utils/chain_type.hpp"
#include "gene_chooser/abstract_gene_chooser.hpp"
#include "nucleotides_remover/abstract_nucleotides_remover.hpp"
#include "p_nucleotides_creator/abstract_nucleotides_creator.hpp"
#include "n_nucleotides_inserter/abstract_n_nucleotides_inserter.hpp"
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

    const double prob_cleavage_v;
    const double prob_cleavage_j;

    const AbstractVDJGeneChooserPtr gene_chooser_p;
    const AbstractNucleotidesRemoverPtr nucl_remover_p;
    const AbstractPNucleotidesCreatorPtr nucl_creator_p;
    const AbstractNNucleotidesInserterPtr nucl_inserter_p;

    const cdr_labeler::DbCDRLabeling v_cdr_db;
    const cdr_labeler::DbCDRLabeling j_cdr_db;

public:
    AbstractMetarootCreator(nullptr_t, nullptr_t,
                            const double prob_cleavage_v,
                            const double prob_cleavage_j,
                            AbstractVDJGeneChooserPtr&& gene_chooser_p_,
                            AbstractNucleotidesRemoverPtr&& nucl_remover_p_,
                            AbstractPNucleotidesCreatorPtr&& nucl_creator_p_,
                            AbstractNNucleotidesInserterPtr&& nucl_inserter_p_,
                            const cdr_labeler::CDRLabelerConfig::CDRsParams &cdr_params_) = delete;

    AbstractMetarootCreator(germline_utils::CustomGeneDatabase *v_db_p,
                            germline_utils::CustomGeneDatabase *j_db_p,
                            const double prob_cleavage_v,
                            const double prob_cleavage_j,
                            AbstractVDJGeneChooserPtr&& gene_chooser_p_,
                            AbstractNucleotidesRemoverPtr&& nucl_remover_p_,
                            AbstractPNucleotidesCreatorPtr&& nucl_creator_p_,
                            AbstractNNucleotidesInserterPtr&& nucl_inserter_p_,
                            const cdr_labeler::CDRLabelerConfig::CDRsParams &cdr_params_) :
        v_db_p(v_db_p),
        j_db_p(j_db_p),
        prob_cleavage_v(prob_cleavage_v),
        prob_cleavage_j(prob_cleavage_j),
        gene_chooser_p(std::move(gene_chooser_p_)),
        nucl_remover_p(std::move(nucl_remover_p_)),
        nucl_creator_p(std::move(nucl_creator_p_)),
        nucl_inserter_p(std::move(nucl_inserter_p_)),
        v_cdr_db(cdr_labeler::GermlineDbLabeler(*v_db_p, cdr_params_).ComputeLabeling()),
        j_cdr_db(cdr_labeler::GermlineDbLabeler(*j_db_p, cdr_params_).ComputeLabeling())
    {
        VERIFY(v_db_p != nullptr);
        VERIFY(j_db_p != nullptr);
        VERIFY(v_db_p->size() > 0);
        VERIFY(j_db_p->size() > 0);
        VERIFY(gene_chooser_p != nullptr);
        VERIFY(nucl_remover_p != nullptr);
        VERIFY(nucl_creator_p != nullptr);
        VERIFY(nucl_inserter_p != nullptr);
        VERIFY(prob_cleavage_v >= 0 and prob_cleavage_v <= 1);
        VERIFY(prob_cleavage_j >= 0 and prob_cleavage_j <= 1);
    }

    virtual ~AbstractMetarootCreator() { }
};

class VJMetarootCreator : public AbstractMetarootCreator {
public:
    VJMetarootCreator(nullptr_t, nullptr_t,
                      const double prob_cleavage_v,
                      const double prob_cleavage_j,
                      AbstractVDJGeneChooserPtr&& gene_chooser_p,
                      AbstractNucleotidesRemoverPtr&& nucl_remover_p,
                      AbstractPNucleotidesCreatorPtr&& nucl_creator_p,
                      AbstractNNucleotidesInserterPtr&& nucl_inserter_p,
                      const cdr_labeler::CDRLabelerConfig::CDRsParams &cdr_params_) = delete;

    VJMetarootCreator(germline_utils::CustomGeneDatabase *v_db_p,
                      germline_utils::CustomGeneDatabase *j_db_p,
                      const double prob_cleavage_v,
                      const double prob_cleavage_j,
                      AbstractVDJGeneChooserPtr&& gene_chooser_p,
                      AbstractNucleotidesRemoverPtr&& nucl_remover_p,
                      AbstractPNucleotidesCreatorPtr&& nucl_creator_p,
                      AbstractNNucleotidesInserterPtr&& nucl_inserter_p,
                      const cdr_labeler::CDRLabelerConfig::CDRsParams &cdr_params_) :
        AbstractMetarootCreator(v_db_p, j_db_p,
                                prob_cleavage_v, prob_cleavage_j,
                                std::move(gene_chooser_p), std::move(nucl_remover_p),
                                std::move(nucl_creator_p), std::move(nucl_inserter_p),
                                cdr_params_)
    { }

    VJMetarootCreator(germline_utils::CustomGeneDatabase &v_db,
                      germline_utils::CustomGeneDatabase &j_db,
                      const double prob_cleavage_v,
                      const double prob_cleavage_j,
                      AbstractVDJGeneChooserPtr&& gene_chooser_p,
                      AbstractNucleotidesRemoverPtr&& nucl_remover_p,
                      AbstractPNucleotidesCreatorPtr&& nucl_creator_p,
                      AbstractNNucleotidesInserterPtr&& nucl_inserter_p,
                      const cdr_labeler::CDRLabelerConfig::CDRsParams &cdr_params_) :
        VJMetarootCreator(&v_db, &j_db,
                          prob_cleavage_v, prob_cleavage_j,
                          std::move(gene_chooser_p), std::move(nucl_remover_p),
                          std::move(nucl_creator_p), std::move(nucl_inserter_p),
                          cdr_params_)
    { }


    VJMetarootCreator(const std::vector<germline_utils::CustomGeneDatabase *>& db,
                      const std::vector<double> prob_cleavages,
                      AbstractVDJGeneChooserPtr&& gene_chooser_p,
                      AbstractNucleotidesRemoverPtr&& nucl_remover_p,
                      AbstractPNucleotidesCreatorPtr&& nucl_creator_p,
                      AbstractNNucleotidesInserterPtr&& nucl_inserter_p,
                      const cdr_labeler::CDRLabelerConfig::CDRsParams &cdr_params_) :
        VJMetarootCreator(db.at(0), db.at(1),
                          std::move(prob_cleavages.at(0)), std::move(prob_cleavages.at(1)),
                          std::move(gene_chooser_p), std::move(nucl_remover_p),
                          std::move(nucl_creator_p), std::move(nucl_inserter_p),
                          cdr_params_)
    {
        VERIFY(db.size() == 2);
        VERIFY(prob_cleavages.size() == 2);
    }

    // TODO virtual (need to return ptr)
    VJMetaRoot CreateRoot() const;
};

class VDJMetarootCreator : public AbstractMetarootCreator {
private:
    germline_utils::CustomGeneDatabase * d_db_p;

    const double prob_cleavage_d_left = 0.5;
    const double prob_cleavage_d_right = 0.5;

public:
    VDJMetarootCreator(nullptr_t, nullptr_t, nullptr_t,
                       const double prob_cleavage_v,
                       const double prob_cleavage_d_left,
                       const double prob_cleavage_d_right,
                       const double prob_cleavage_j,
                       AbstractVDJGeneChooserPtr&& gene_chooser_p_,
                       AbstractNucleotidesRemoverPtr&& nucl_remover_p_,
                       AbstractPNucleotidesCreatorPtr&& nucl_creator_p_,
                       AbstractNNucleotidesInserterPtr&& nucl_inserter_p_,
                       const cdr_labeler::CDRLabelerConfig::CDRsParams &cdr_params_) = delete;

    VDJMetarootCreator(germline_utils::CustomGeneDatabase * v_db_p,
                       germline_utils::CustomGeneDatabase * d_db_p,
                       germline_utils::CustomGeneDatabase * j_db_p,
                       const double prob_cleavage_v,
                       const double prob_cleavage_d_left,
                       const double prob_cleavage_d_right,
                       const double prob_cleavage_j,
                       AbstractVDJGeneChooserPtr&& gene_chooser_p_,
                       AbstractNucleotidesRemoverPtr&& nucl_remover_p_,
                       AbstractPNucleotidesCreatorPtr&& nucl_creator_p_,
                       AbstractNNucleotidesInserterPtr&& nucl_inserter_p_,
                       const cdr_labeler::CDRLabelerConfig::CDRsParams &cdr_params_) :
            AbstractMetarootCreator(v_db_p, j_db_p,
                                    prob_cleavage_v, prob_cleavage_j,
                                    std::move(gene_chooser_p_), std::move(nucl_remover_p_),
                                    std::move(nucl_creator_p_), std::move(nucl_inserter_p_),
                                    cdr_params_),
            d_db_p(d_db_p),
            prob_cleavage_d_left(prob_cleavage_d_left),
            prob_cleavage_d_right(prob_cleavage_d_right)
    {
        VERIFY(d_db_p != nullptr);
        VERIFY(d_db_p->size() > 0);
        VERIFY(prob_cleavage_d_left >= 0 and prob_cleavage_d_left <= 1);
        VERIFY(prob_cleavage_d_right >= 0 and prob_cleavage_d_right <= 1);
    }

    VDJMetarootCreator(germline_utils::CustomGeneDatabase &v_db,
                       germline_utils::CustomGeneDatabase &d_db,
                       germline_utils::CustomGeneDatabase &j_db,
                       const double prob_cleavage_v,
                       const double prob_cleavage_d_left,
                       const double prob_cleavage_d_right,
                       const double prob_cleavage_j,
                       AbstractVDJGeneChooserPtr&& gene_chooser_p,
                       AbstractNucleotidesRemoverPtr&& nucl_remover_p,
                       AbstractPNucleotidesCreatorPtr&& nucl_creator_p,
                       AbstractNNucleotidesInserterPtr&& nucl_inserter_p,
                       const cdr_labeler::CDRLabelerConfig::CDRsParams &cdr_params_) :
            VDJMetarootCreator(&v_db, &d_db, &j_db,
                               prob_cleavage_v, prob_cleavage_j,
                               prob_cleavage_d_left, prob_cleavage_d_right,
                               std::move(gene_chooser_p), std::move(nucl_remover_p),
                               std::move(nucl_creator_p), std::move(nucl_inserter_p),
                               cdr_params_)
    { }


    VDJMetarootCreator(const std::vector<germline_utils::CustomGeneDatabase *>& db,
                       std::vector<double> &&prob_cleavages,
                       AbstractVDJGeneChooserPtr&& gene_chooser_p,
                       AbstractNucleotidesRemoverPtr&& nucl_remover_p,
                       AbstractPNucleotidesCreatorPtr&& nucl_creator_p,
                       AbstractNNucleotidesInserterPtr&& nucl_inserter_p,
                       const cdr_labeler::CDRLabelerConfig::CDRsParams &cdr_params_) :
        VDJMetarootCreator(db.at(0), db.at(1), db.at(2),
                           std::move(prob_cleavages.at(0)), std::move(prob_cleavages.at(1)),
                           std::move(prob_cleavages.at(2)), std::move(prob_cleavages.at(3)),
                           std::move(gene_chooser_p), std::move(nucl_remover_p),
                           std::move(nucl_creator_p), std::move(nucl_inserter_p),
                           cdr_params_)
    {
        VERIFY(db.size() == 3);
        VERIFY(prob_cleavages.size() == 4);
    }

    // TODO virtual (need to return ptr)
    VDJMetaRoot CreateRoot() const;
};

} // End namespace ig_simulator

//
// Created by Andrew Bzikadze on 3/16/17.
//

#pragma once

#include <tuple>
#include <memory>

#include "germline_utils/germline_databases/custom_gene_database.hpp"
#include "ig_simulator_utils.hpp"

namespace ig_simulator {

using VDJ_GenesIndexTuple = std::tuple<size_t, size_t, size_t>;

class AbstractVDJGeneChooser {
protected:
    const germline_utils::CustomGeneDatabase *v_db_p_;
    const germline_utils::CustomGeneDatabase *d_db_p_;
    const germline_utils::CustomGeneDatabase *j_db_p_;

    // This variable defines whether D segment is generated
    // d_dp_p_ MUST be nullptr if is_vdj == false
    bool is_vdj;

public:
    explicit AbstractVDJGeneChooser(const std::vector<germline_utils::CustomGeneDatabase>& db);

    virtual VDJ_GenesIndexTuple ChooseGenes() const = 0;

    /**
     * This method suggests whether D segment is generated.
     * If `false` then second component of `VDJ_GenesIndexTuple`
     * returned by `ChooseGenes()` will be size_t(-1).
     */
    bool IsVDJ() const;

    AbstractVDJGeneChooser() = delete;
    AbstractVDJGeneChooser(const AbstractVDJGeneChooser&) = delete;
    AbstractVDJGeneChooser(AbstractVDJGeneChooser&&) = delete;
    AbstractVDJGeneChooser& operator=(const AbstractVDJGeneChooser&) = delete;
    AbstractVDJGeneChooser& operator=(AbstractVDJGeneChooser&&) = delete;

    virtual ~AbstractVDJGeneChooser() { }
};

using AbstractVDJGeneChooserCPtr = std::unique_ptr<const AbstractVDJGeneChooser>;

} // End namespace ig_simulator

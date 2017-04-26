//
// Created by Andrew Bzikadze on 3/16/17.
//

#include "abstract_gene_chooser.hpp"

namespace ig_simulator {

AbstractVDJGeneChooser::AbstractVDJGeneChooser(const std::vector<germline_utils::CustomGeneDatabase>& db):
    v_db_p_(&db.front()),
    d_db_p_(nullptr),
    j_db_p_(&db.back()),
    is_vdj(false)
{
    VERIFY(db.size() >= 2 and db.size() <= 3);

    if (db.size() == 3) {
        d_db_p_ = &db[1];
        is_vdj = true;
    }
}

inline bool AbstractVDJGeneChooser::IsVDJ() const {
    // if (not is_vdj) { VERIFY(d_db_p_ != nullptr); }
    return is_vdj;
}

} // End namespace ig_simulator

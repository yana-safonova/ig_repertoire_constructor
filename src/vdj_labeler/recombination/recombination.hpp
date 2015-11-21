#pragma once

#include "cleaved_gene.hpp"
#include "nongenomic_insertion.hpp"

class HCRecombination {
    // here should be reference to read
    const CleavedIgGeneAlignment& v_gene_;
    const CleavedIgGeneAlignment& d_gene_;
    const CleavedIgGeneAlignment& j_gene_;
    NongenomicInsertion vd_insertion_;
    NongenomicInsertion dj_insertion_;

public:
    HCRecombination(const CleavedIgGeneAlignment& v_gene,
                    const CleavedIgGeneAlignment& d_gene,
                    const CleavedIgGeneAlignment& j_gene,
                    NongenomicInsertion vd_insertion,
                    NongenomicInsertion dj_insertion) :
            v_gene_(v_gene),
            d_gene_(d_gene),
            j_gene_(j_gene),
            vd_insertion_(vd_insertion),
            dj_insertion_(dj_insertion) { }

    const CleavedIgGeneAlignment& V() const { return v_gene_; }

    const CleavedIgGeneAlignment& D() const { return d_gene_; }

    const CleavedIgGeneAlignment& J() const { return j_gene_; }

    NongenomicInsertion VDInsertion() const { return vd_insertion_; }

    NongenomicInsertion DJInsertion() const { return dj_insertion_; }
};
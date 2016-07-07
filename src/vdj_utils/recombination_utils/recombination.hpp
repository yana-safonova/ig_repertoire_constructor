#pragma once

#include "cleaved_gene.hpp"
#include "nongenomic_insertion.hpp"
#include <read_archive.hpp>

namespace recombination_utils {

class HCRecombination {
    core::ReadPtr read_ptr_;
    CleavedIgGeneAlignment v_gene_;
    CleavedIgGeneAlignment d_gene_;
    CleavedIgGeneAlignment j_gene_;
    NongenomicInsertion vd_insertion_;
    NongenomicInsertion dj_insertion_;

public:
    HCRecombination(core::ReadPtr read_ptr,
                    CleavedIgGeneAlignment v_gene,
                    CleavedIgGeneAlignment d_gene,
                    CleavedIgGeneAlignment j_gene,
                    NongenomicInsertion vd_insertion,
                    NongenomicInsertion dj_insertion) :
        read_ptr_(read_ptr),
        v_gene_(v_gene),
        d_gene_(d_gene),
        j_gene_(j_gene),
        vd_insertion_(vd_insertion),
        dj_insertion_(dj_insertion) { }

    HCRecombination(const HCRecombination &obj) :
        v_gene_(obj.v_gene_),
        d_gene_(obj.d_gene_),
        j_gene_(obj.j_gene_) {
        read_ptr_ = obj.read_ptr_;
        vd_insertion_ = obj.vd_insertion_;
        dj_insertion_ = obj.dj_insertion_;
    }

    core::ReadPtr Read() const { return read_ptr_; }

    const CleavedIgGeneAlignment &V() const { return v_gene_; }

    const CleavedIgGeneAlignment &D() const { return d_gene_; }

    const CleavedIgGeneAlignment &J() const { return j_gene_; }

    NongenomicInsertion VDInsertion() const { return vd_insertion_; }

    NongenomicInsertion DJInsertion() const { return dj_insertion_; }

    size_t SHMsNumber() const { return v_gene_.SHMsNumber() + d_gene_.SHMsNumber() + j_gene_.SHMsNumber(); }

    bool Valid() const;
};

std::ostream &operator<<(std::ostream &out, const HCRecombination &hc_recombination);

} // End namespace recombination_utils
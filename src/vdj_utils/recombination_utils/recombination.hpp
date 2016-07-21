#pragma once

#include "cleaved_gene.hpp"
#include "nongenomic_insertion.hpp"
#include <read_archive.hpp>

namespace recombination_utils {

class HCRecombination {
    const core::Read* read_ptr_;
    CleavedIgGeneAlignment v_gene_;
    CleavedIgGeneAlignment d_gene_;
    CleavedIgGeneAlignment j_gene_;
    NongenomicInsertion vd_insertion_;
    NongenomicInsertion dj_insertion_;

public:
    HCRecombination() : read_ptr_(nullptr),
                        v_gene_(), d_gene_(), j_gene_(),
                        vd_insertion_(), dj_insertion_()
    { }

    HCRecombination(const core::Read* read_ptr,
                    CleavedIgGeneAlignment v_gene,
                    CleavedIgGeneAlignment d_gene,
                    CleavedIgGeneAlignment j_gene,
                    NongenomicInsertion vd_insertion,
                    NongenomicInsertion dj_insertion) :
        read_ptr_(read_ptr),
        v_gene_(std::move(v_gene)),
        d_gene_(std::move(d_gene)),
        j_gene_(std::move(j_gene)),
        vd_insertion_(std::move(vd_insertion)),
        dj_insertion_(std::move(dj_insertion))
    { }

    HCRecombination(const core::Read& read,
                    CleavedIgGeneAlignment v_gene,
                    CleavedIgGeneAlignment d_gene,
                    CleavedIgGeneAlignment j_gene,
                    NongenomicInsertion vd_insertion,
                    NongenomicInsertion dj_insertion) :
        HCRecombination(&read,
                        std::move(v_gene), std::move(d_gene), std::move(j_gene),
                        std::move(vd_insertion), std::move(dj_insertion))
    { }

    HCRecombination(const HCRecombination&) = default;
    HCRecombination(HCRecombination&&) = default;
    HCRecombination& operator=(const HCRecombination&) = default;
    HCRecombination& operator=(HCRecombination&&) = default;

    const core::Read* ReadPtr() const { return read_ptr_; }
    const core::Read& Read() const {
        assert(read_ptr_ != nullptr);
        return *read_ptr_;
    }

    const CleavedIgGeneAlignment &V() const { return v_gene_; }
    const CleavedIgGeneAlignment &D() const { return d_gene_; }
    const CleavedIgGeneAlignment &J() const { return j_gene_; }

    const NongenomicInsertion &VDInsertion() const { return vd_insertion_; }
    const NongenomicInsertion &DJInsertion() const { return dj_insertion_; }

    size_t SHMsNumber() const { return v_gene_.SHMsNumber() + d_gene_.SHMsNumber() + j_gene_.SHMsNumber(); }

    // recombination is not valid if v gene event overlaps j gene event or
    // j gene event overlaps d gene event
    // it result in invalid vd/dj insertion
    bool Valid() const {
        return vd_insertion_.Valid() and dj_insertion_.Valid();
    }

    size_t ReadId() const {
        assert(read_ptr_ != nullptr);
        return read_ptr_-> id;
    }
};

std::ostream &operator<<(std::ostream &out, const HCRecombination &hc_recombination);

} // End namespace recombination_utils
//
// Created by Andrew Bzikadze on 5/21/16.
//

#pragma once

#include <memory>
#include "gene_alignment/gene_alignment.hpp"

namespace ns_abstract_alignment_cropper {
class AbstractAlignmentCropper {
public:
    virtual void crop(ns_gene_alignment::ReadGermlineAlignment&) const = 0;
    virtual ~AbstractAlignmentCropper() { }
};
using AbstractAlignmentCropperPtr = std::shared_ptr<AbstractAlignmentCropper>;
}

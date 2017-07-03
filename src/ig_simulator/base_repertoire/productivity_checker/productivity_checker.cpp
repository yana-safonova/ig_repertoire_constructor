//
// Created by Andrew Bzikadze on 3/27/17.
//

#include "productivity_checker.hpp"

namespace ig_simulator {

bool ProductivityChecker::IsProductive(const AbstractMetaroot& root) const {
    if (root.CDRLabeling().Empty())
        return false;
    core::Read read("", root.Sequence(), 0);
    auto aa = aa_calculator->ComputeAminoAcidAnnotation(read, root.CDRLabeling());
    return not aa.HasStopCodon() and aa.InFrame();
}

} // End namespace ig_simulator
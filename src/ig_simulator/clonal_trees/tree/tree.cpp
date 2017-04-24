//
// Created by Andrew Bzikadze on 4/9/17.
//

#include "tree.hpp"
#include "verify.hpp"
#include "annotation_utils/aa_annotation/aa_calculator.hpp"

namespace ig_simulator {

std::ostream& operator<<(std::ostream& out, const Tree& tree) {
    VERIFY(tree.nodes.size() >= 1);

    for (size_t i = 1; i < tree.nodes.size(); ++i) {
        const auto& node = tree.nodes[i];
        const auto& shms = node.SHMs();
        out << node.ParentInd() << " -> " << i << ' ';
        out << '(' << (node.IsIncluded() ? "included" : "excluded") << ')' << ' ';
        out << " shms: ";
        for(const auto& shm : shms) {
            out << "(at " << std::get<0>(shm) <<
                   ", from " << std::get<1>(shm) <<
                   ", to " << std::get<2>(shm) << ')' << ' ';
        }
        out << '\n';

        // annotation_utils::SimpleAACalculator aa_calculator;
        // core::Read read("", tree.Sequences()[node.ParentInd()], 0);
        // VERIFY(not aa_calculator.ComputeAminoAcidAnnotation(read, tree.Metaroot()->CDRLabeling()).HasStopCodon());
    }
    return out;
}

} // End namespace ig_simulator

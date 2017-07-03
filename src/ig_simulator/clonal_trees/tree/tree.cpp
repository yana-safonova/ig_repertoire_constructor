//
// Created by Andrew Bzikadze on 4/9/17.
//

#include "tree.hpp"
#include "verify.hpp"
#include "annotation_utils/aa_annotation/aa_calculator.hpp"

namespace ig_simulator {

std::ostream& operator<<(std::ostream& out, const Tree& tree) {
    VERIFY(tree.nodes.size() >= 1);

    out << "digraph G {\n";
    out << '\t' << 0 << " [shape = "     << (tree.Metaroot()->IsProductive() ? "circle" : "box") << "," <<
                          "fillcolor = " << (tree.nodes.front().IsIncluded() ? "cyan"   : "magenta") << "," <<
                          "style = filled,size=1]; // " <<
        '(' << (tree.nodes.front().IsIncluded()   ? "included" : "excluded") << ')' << ' ' <<
        '(' << (tree.nodes.front().IsProductive() ? "productive" : "non-productive") << ')' << '\n';

    for (size_t i = 1; i < tree.nodes.size(); ++i) {
        const auto& node = tree.nodes[i];
        const auto& shms = node.SHMs();
        out << '\t' << i << " [shape = "     << (node.IsProductive() ? "circle" : "box"    ) << "," <<
                              "fillcolor = " << (node.IsIncluded()   ? "cyan"   : "magenta") << "," <<
                              "style = filled,size=1]; // " <<
            '(' << (node.IsIncluded() ? "included" : "excluded") << ')' << ' ' <<
            '(' << (node.IsProductive() ? "productive" : "non-productive") << ')' << '\n';
        out << '\t' << node.ParentInd() << " -> " << i << "[minlen = " << std::to_string(shms.size()) << "]; // ";
        out << "total shms: " << shms.size() << " ";
        for(const auto& shm : shms) {
            out << "(at " << std::get<0>(shm) <<
                   " from " << std::get<1>(shm) <<
                   " to " << std::get<2>(shm) << ')' << ' ';
        }
        out << '\n';

        // annotation_utils::SimpleAACalculator aa_calculator;
        // core::Read read("", tree.Sequences()[node.ParentInd()], 0);
        // VERIFY(not aa_calculator.ComputeAminoAcidAnnotation(read, tree.Metaroot()->CDRLabeling()).HasStopCodon());
    }
    out << "}\n";
    return out;
}

} // End namespace ig_simulator

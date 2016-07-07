#include "cleaved_gene.hpp"

namespace recombination_utils {

std::ostream& operator<<(std::ostream &out, const CleavedIgGeneAlignment& cleaved_gene) {
    out << cleaved_gene.GeneAlignment() << std::endl;
    out << "Left cleavage len: " << cleaved_gene.LeftCleavageLength() <<
            ", right cleavage len: " << cleaved_gene.RightCleavageLength() << std::endl;
    out << "Start position: " << cleaved_gene.StartReadPosition() <<
            ", end position: " << cleaved_gene.EndReadPosition() << std::endl;
    out << "# SHMs: " << cleaved_gene.SHMsNumber() << std::endl;
    return out;
}

} // End namespace recombination_utils
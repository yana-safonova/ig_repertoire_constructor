#include "recombination.hpp"

namespace recombination_utils {


std::ostream &operator<<(std::ostream &out, const HCRecombination &hc_recombination) {
    out << "HC recombination, # SHMs: " << hc_recombination.SHMsNumber() << ", valid: " <<
        hc_recombination.Valid() << std::endl;
    out << "V event: " << std::endl;
    out << hc_recombination.V() << std::endl;
    out << "D event: " << std::endl;
    out << hc_recombination.D() << std::endl;
    out << "J event: " << std::endl;
    out << hc_recombination.J() << std::endl;
    out << "VD insertion: " << hc_recombination.VDInsertion() << std::endl;
    out << "DJ insertion: " << hc_recombination.DJInsertion() << std::endl;
    return out;
}

} // End namespace recombination_utils
#include "lymphocyte_type.hpp"

namespace germline_utils {

    std::ostream &operator<<(std::ostream &out, const LymphocyteType &lymphocyte_type) {
        if (lymphocyte_type == LymphocyteType::BLymphocyte)
            out << "B-cell";
        else if (lymphocyte_type == LymphocyteType::TLymphocyte)
            out << "T-cell";
        else
            out << "unknown";
        return out;
    }

}
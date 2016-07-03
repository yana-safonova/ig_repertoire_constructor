#pragma once

#include <ostream>

namespace germline_utils {

    enum LymphocyteType {
        UnknownLymphocyte, BLymphocyte, TLymphocyte
    };

    std::ostream &operator<<(std::ostream &out, const LymphocyteType &lymphocyte_type);

}
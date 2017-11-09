#pragma once

#include <boost/optional/optional.hpp>
#include "../ig_tools/utils/string_tools.hpp"

struct CloneInfo {
    std::size_t id;
    std::size_t size;

    CloneInfo(std::size_t id, std::size_t size) : id(id), size(size) {}

    static boost::optional<CloneInfo> TryParse(const std::string& name);
};

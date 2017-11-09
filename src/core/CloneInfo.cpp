#include "CloneInfo.hpp"

boost::optional<CloneInfo> CloneInfo::TryParse(const std::string& name) {
    const std::vector<std::string> parts = split(name, "___");
    if (parts.size() != 4) return {};
    const auto id = try_string_to_number<std::size_t>(parts[1]);
    const auto size = try_string_to_number<std::size_t>(parts[3]);
    if (!id || !size) return {};
    return {CloneInfo(*id, *size)};
}

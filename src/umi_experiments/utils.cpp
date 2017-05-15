#include "utils.hpp"

void create_console_logger(logging::level log_level) {
    using namespace logging;
    logger *lg = create_logger("", log_level);
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}
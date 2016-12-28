#include "utils.hpp"
#include <logger/log_writers.hpp>

using bformat = boost::format;


std::string join_cmd_line(size_t argc, char **argv) {
    std::string result = argv[0];
    for (size_t i = 1; i < argc; ++i) {
        result += " ";
        result += argv[i];
    }

    return result;
}


void create_console_logger(std::string log_props_file) {
    using namespace logging;

    logger *lg = create_logger(log_props_file);
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}


std::string running_time_format(const perf_counter &pc) {
    unsigned ms = (unsigned)pc.time_ms();
    unsigned secs = (ms / 1000) % 60;
    unsigned mins = (ms / 1000 / 60) % 60;
    unsigned hours = (ms / 1000 / 60 / 60);
    bformat bf("%u hours %u minutes %u seconds");
    bf % hours % mins % secs;
    return bf.str();
}

#pragma once

#include <perfcounter.hpp>

std::string join_cmd_line(size_t argc, char **argv);

void create_console_logger(std::string log_props_file = "");

std::string running_time_format(const perf_counter &pc);

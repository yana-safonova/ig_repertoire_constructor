#pragma once

#include <logger/log_writers.hpp>
#include <seqan/seq_io.h>

void create_console_logger();

size_t get_sw_dist(const seqan::Dna5String& first, const seqan::Dna5String& second);

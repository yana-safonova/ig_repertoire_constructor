#pragma once

#include <seqan/seq_io.h>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem.hpp>

bool create_new_directory(boost::filesystem::path path);

void write_seqan_records(boost::filesystem::path path, std::vector<seqan::CharString> ids, std::vector<seqan::Dna5String> reads);
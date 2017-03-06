#pragma once

#include <unordered_map>
#include <seqan/seq_io.h>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem.hpp>

bool create_new_directory(boost::filesystem::path path);

void write_seqan_records(boost::filesystem::path path, std::vector<seqan::CharString> ids, std::vector<seqan::Dna5String> reads);

void read_seqan_records(const std::string& input_file_name, std::vector<seqan::CharString>& ids, std::vector<seqan::Dna5String>& reads);

std::unordered_map<seqan::CharString, size_t> read_rcm_file(std::string file_path);
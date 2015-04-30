#pragma once

#include <vector>
#include "spliced_read.hpp"

namespace ig_repertoire_constructor {

class SplicerResult {
public:
    std::vector <size_t> spliced_read_indices;
	std::vector <SplicedRead> spliced_reads;

	std::vector <size_t> bad_read_numbers;
};

}

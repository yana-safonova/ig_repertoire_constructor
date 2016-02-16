#pragma once

#include "abpair_config.hpp"

class AbPairLauncher {

    std::vector<std::string> ReadInputFnames(std::string input_sequences);

public:
    void Run(const abpair_config::io_config &io);
};
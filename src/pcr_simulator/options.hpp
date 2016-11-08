#pragma once

struct Options {
    struct SimulationOptions {
        size_t cycles_count;
        double error_prob_first;
        double error_prob_last;
        double amplification_rate;
        double chimeras_rate;
        // 1 for barcode going with the left half, 2 for the right half, 3 for random choice
        size_t barcode_position;
        size_t barcode_length;
    };

    std::string repertoire_file_path;
    std::string output_dir_path;
    size_t output_estimation_limit;
    SimulationOptions simulation_options;
};
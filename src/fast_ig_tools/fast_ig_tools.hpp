#pragma once

#include <string>
#include <path_helper.hpp>
#include <perfcounter.hpp>

#include <logger/log_writers.hpp>
#include <logger/logger.hpp>
#include <segfault_handler.hpp>

#include <boost/format.hpp>
using bformat = boost::format;

using Graph = std::vector<std::vector<std::pair<size_t, int>>>;

size_t numEdges(const Graph &graph,
                bool undirected = true);

void write_metis_graph(const Graph &graph,
                       const std::string &filename,
                       bool undirected = true);

void write_metis_graph(const Graph &graph,
                       const std::vector<size_t> &weights,
                       const std::string &filename,
                       bool undirected = true);

std::vector<size_t> optimal_coverage(const std::vector<size_t> &multiplicities,
                                     size_t K, size_t n = 3);

// vim: ts=4:sw=4

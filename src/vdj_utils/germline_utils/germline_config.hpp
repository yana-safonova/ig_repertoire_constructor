//
// Created by Andrew Bzikadze on 3/16/17.
//

#pragma once
#include "io/library.hpp"
#include "config_singl.hpp"
#include <boost/property_tree/ptree_fwd.hpp>

namespace germline_utils {

struct GermlineInput {
    std::string ig_dir;
    std::string tcr_dir;
    std::string germline_filenames_config;
};

struct GermlineParams {
    std::string germline_dir;
    std::string organism;
    std::string loci;
    bool pseudogenes;
};


void load(GermlineInput &gi, boost::property_tree::ptree const &pt, bool);
void load(GermlineParams &gp, boost::property_tree::ptree const &pt, bool);

} // End namespace germline_utils

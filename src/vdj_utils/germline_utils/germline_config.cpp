//
// Created by Andrew Bzikadze on 3/16/17.
//

#include "germline_config.hpp"
#include <logger/logger.hpp>
#include <config_common.hpp>
#include <boost/property_tree/ptree_fwd.hpp>

namespace germline_utils {

void load(GermlineInput &gi, boost::property_tree::ptree const &pt, bool) {
    using config_common::load;
    load(gi.germline_filenames_config, pt, "germline_filenames_config");
    load(gi.ig_dir, pt, "ig_dir");
    load(gi.tcr_dir, pt, "tcr_dir");
}

void load(GermlineParams &gp, boost::property_tree::ptree const &pt, bool) {
    using config_common::load;
    load(gp.germline_dir, pt, "germline_dir");
    load(gp.loci, pt, "loci");
    load(gp.organism, pt, "organism");
    load(gp.pseudogenes, pt, "pseudogenes");
}


} // End namespace germline_utils

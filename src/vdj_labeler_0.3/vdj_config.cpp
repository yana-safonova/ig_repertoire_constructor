#include "vdj_config.hpp"
#include "../include/config_common.hpp"

namespace vdj_labeler {

void load(VDJLabelerConfig::IOParams::InputParams::GermlineGenes::IGH_Genes &igh,
          boost::property_tree::ptree const &pt, bool) {
    using config_common::load;
    load(igh.variable_genes, pt, "variable_genes");
    load(igh.diversity_genes, pt, "diversity_genes");
    load(igh.join_genes, pt, "join_genes");
}

void load(VDJLabelerConfig::IOParams::InputParams::GermlineGenes &germlines,
          boost::property_tree::ptree const &pt, bool) {
    using config_common::load;
    load(germlines.igh_genes, pt, "igh");
}

void load(VDJLabelerConfig::IOParams::InputParams &input, boost::property_tree::ptree const &pt, bool) {
    using config_common::load;
    load(input.input_sequences, pt, "input_sequences");
    load(input.vj_finder_config, pt, "vj_finder_config");
    load(input.germlines, pt, "germlines");
}

void load(VDJLabelerConfig::IOParams::OutputParams &output, boost::property_tree::ptree const &pt, bool) {
    using config_common::load;
    load(output.log_filename, pt, "log_filename");
    load(output.output_dir, pt, "output_dir");
}

void load(VDJLabelerConfig::IOParams &iop, boost::property_tree::ptree const &pt, bool) {
    using config_common::load;
    load(iop.input_params, pt, "input_params");
    load(iop.output_params, pt, "output_params");
}

void load(VDJLabelerConfig::RunParams &rp, boost::property_tree::ptree const &pt, bool) {
    using config_common::load;
    load(rp.max_memory, pt, "max_memory");
    load(rp.threads_count, pt, "thread_count");
}

void VDJLabelerConfig::load(std::string const &filename) {
    boost::property_tree::ptree pt;
    boost::property_tree::read_info(filename, pt);
    using config_common::load;
    load(io_params, pt, "io_params");
    load(run_params, pt, "run_params");
    vj_finder::load(vj_finder_config, io_params.input_params.vj_finder_config);
}

} // End namespace vdj_labeler
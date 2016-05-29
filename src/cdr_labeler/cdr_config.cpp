#include "cdr_config.hpp"
#include <config_common.hpp>

namespace cdr_labeler {

    void load(CDRLabelerConfig::InputParams &ip, boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(ip.input_reads, pt, "input_reads");
        load(ip.vj_finder_config, pt, "vj_finder_config");
    }

    void load(CDRLabelerConfig::OutputParams &op, boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(op.output_dir, pt, "output_dir");
    }

    void load(CDRLabelerConfig::RunParams &rp, boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(rp.num_threads, pt, "num_threads");
    }

    CDRLabelerConfig::CDRsParams::AnnotatedSearchParams::DomainSystem convert_string_to_domain_system(std::string str) {
        if(str == "imgt")
            return CDRLabelerConfig::CDRsParams::AnnotatedSearchParams::DomainSystem::IMGT;
        if(str == "kabat")
            return CDRLabelerConfig::CDRsParams::AnnotatedSearchParams::DomainSystem::Kabat;
        VERIFY_MSG(false, "Domain system " << str << " is unknown!");
        return CDRLabelerConfig::CDRsParams::AnnotatedSearchParams::DomainSystem::UnknownDomain;
    }

    void load(CDRLabelerConfig::CDRsParams::AnnotatedSearchParams &ap, boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        std::string domain_system_str;
        load(domain_system_str, pt, "domain_system");
        ap.domain_system = convert_string_to_domain_system(domain_system_str);
        load(ap.imgt_j_annotation, pt, "imgt_j_annotation");
        load(ap.kabat_j_annotation, pt, "kabat_j_annotation");
        load(ap.imgt_v_annotation, pt, "imgt_v_annotation");
        load(ap.kabat_v_annotation, pt, "kabat_v_annotation");
    }

    void load(CDRLabelerConfig::CDRsParams::HCDR1Params &cdr1_p, boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(cdr1_p.max_length, pt, "max_length");
        load(cdr1_p.min_length, pt, "min_length");
        load(cdr1_p.residues_after, pt, "residues_after");
        load(cdr1_p.residues_before, pt, "residues_before");
        load(cdr1_p.start_pos, pt, "start_pos");
        load(cdr1_p.start_shift, pt, "start_shift");
    }

    void load(CDRLabelerConfig::CDRsParams::HCDR2Params &cdr2_p, boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(cdr2_p.max_length, pt, "max_length");
        load(cdr2_p.min_length, pt, "min_length");
        load(cdr2_p.residues_after, pt, "residues_after");
        load(cdr2_p.residues_before, pt, "residues_before");
        load(cdr2_p.distance_from_cdr1_end, pt, "distance_from_cdr1_end");
        load(cdr2_p.distance_shift, pt, "distance_shift");
    }

    void load(CDRLabelerConfig::CDRsParams::HCDR3Params &cdr3_p, boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(cdr3_p.max_length, pt, "max_length");
        load(cdr3_p.min_length, pt, "min_length");
        load(cdr3_p.residues_after, pt, "residues_after");
        load(cdr3_p.residues_before, pt, "residues_before");
        load(cdr3_p.distance_from_cdr2_end, pt, "distance_from_cdr2_end");
        load(cdr3_p.distance_shift, pt, "distance_shift");
    }

    CDRLabelerConfig::CDRsParams::CDRSearchAlgorithm convert_str_cdr_search_params(std::string str) {
        if(str == "annotated_search")
            return CDRLabelerConfig::CDRsParams::CDRSearchAlgorithm::AnnotatedSearch;
        if(str == "de_novo_search")
            return CDRLabelerConfig::CDRsParams::CDRSearchAlgorithm::DeNovoSearch;
        VERIFY_MSG(false, "Unknown CDR search algorithm: " << str);
        return CDRLabelerConfig::CDRsParams::CDRSearchAlgorithm::UnknownSearchAlgorithm;
    }

    void load(CDRLabelerConfig::CDRsParams &cdrs_p, boost::property_tree::ptree const &pt, bool) {
        using config_common::load;
        load(cdrs_p.annotated_search_params, pt, "annotated_search_params");
        load(cdrs_p.hcdr1_params, pt, "hcdr1_params");
        load(cdrs_p.hcdr2_params, pt, "hcdr2_params");
        load(cdrs_p.hcdr3_params, pt, "hcdr3_params");
        std::string cdr_search_str;
        load(cdr_search_str, pt, "cdr_search_algorithm");
        cdrs_p.cdr_search_algorithm = convert_str_cdr_search_params(cdr_search_str);
    }

    void CDRLabelerConfig::load(std::string config_fname) {
        boost::property_tree::ptree pt;
        boost::property_tree::read_info(config_fname, pt);
        using config_common::load;
        load(input_params, pt, "input_params", true);
        load(output_params, pt, "output_params", true);
        load(run_params, pt, "run_params", true);
        load(cdrs_params, pt, "cdrs_params", true);
        vj_finder::load(vj_finder_config, input_params.vj_finder_config);
        vj_finder_config.algorithm_params.germline_params.pseudogenes = false;
        vj_finder_config.algorithm_params.germline_params.loci = "IGH";
    }
}
#pragma once

#include "vj_finder_config.hpp"

namespace cdr_labeler {
    struct CDRLabelerConfig {
        struct InputParams {
            std::string input_reads;
            std::string vj_finder_config;
            std::string run_hg_constructor;
        };

        struct OutputParams {
            struct FeatureReportParams {
                std::string cdr_details;
                std::string preset;
                std::string columns;
            };

            std::string output_dir;
            std::string shm_details;
            std::string cdr1_fasta;
            std::string cdr2_fasta;
            std::string cdr3_fasta;
            std::string cdr3_compressed_fasta;
            std::string v_alignment_fasta;
            std::string cleaned_reads;
            std::string trash_output;
            FeatureReportParams feature_report_params;
        };

        struct RunParams {
            size_t num_threads;
        };

        struct CDRsParams {
            struct AnnotatedSearchParams {
                enum DomainSystem { Unknown_Domain, IMGT_Domain, Kabat_Domain };

                struct VGeneAnnotation {
                    std::string imgt_v_annotation;
                    std::string kabat_v_annotation;
                    size_t v_gene_line_index;
                    size_t cdr1_start_line_index;
                    size_t cdr1_end_line_index;
                    size_t cdr2_start_line_index;
                    size_t cdr2_end_line_index;
                    size_t fr3_end_index;
                };

                struct JGeneAnnotation {
                    std::string imgt_j_annotation;
                    std::string kabat_j_annotation;
                    size_t j_gene_line_index;
                    size_t cdr3_end_index;
                };

                DomainSystem domain_system;
                VGeneAnnotation v_gene_annotation;
                JGeneAnnotation j_gene_annotation;
            };

            struct HCDR1Params {
                size_t start_pos;
                size_t start_shift;
                size_t min_length;
                size_t max_length;
                std::string residues_before;
                std::string residues_after;
            };

            struct HCDR2Params {
                size_t distance_from_cdr1_end;
                size_t distance_shift;
                size_t min_length;
                size_t max_length;
                std::string residues_before;
                std::string residues_after;
            };

            struct HCDR3Params {
                size_t distance_from_cdr2_end;
                size_t distance_shift;
                size_t min_length;
                size_t max_length;
                std::string residues_before;
                std::string residues_after;
            };

            AnnotatedSearchParams annotated_search_params;
            HCDR1Params hcdr1_params;
            HCDR2Params hcdr2_params;
            HCDR3Params hcdr3_params;

            enum CDRSearchAlgorithm { UnknownCDRSearchAlgorithm, AnnotatedCDRSearch, DeNovoCDRSearch };
            CDRSearchAlgorithm cdr_search_algorithm;
        };

        struct SHMFindingParams {
            struct SHMFilteringParams {
                size_t v_start_max_skipped;
                size_t v_end_max_skipped;
                size_t j_start_max_skipped;
                size_t j_end_max_skipped;
            };

            enum SHMFindingAlgorithm { UnknownSHMSearchAlgorithm, FilteringSHMAlgorithm, CDRFilteringSHMAlgorithm };

            SHMFilteringParams shm_filtering_params;
            SHMFindingAlgorithm shm_finding_algorithm;
        };

        InputParams input_params;
        OutputParams output_params;
        RunParams run_params;
        CDRsParams cdrs_params;
        SHMFindingParams shm_params;

        vj_finder::VJFinderConfig vj_finder_config;

    };

    void load(CDRLabelerConfig& cfg, const std::string& config_fname);

    typedef config_common::config<CDRLabelerConfig> cdrl_cfg;
}

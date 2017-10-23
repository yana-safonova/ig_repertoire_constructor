#include <gtest/gtest.h>
#include <logger/log_writers.hpp>

#include <sstream>

#include <cdr_config.hpp>
#include <germline_utils/germline_db_generator.hpp>
#include <germline_db_labeler.hpp>
#include <vj_parallel_processor.hpp>
#include <read_labeler.hpp>
#include <convert.hpp>
#include <annotation_utils/shm_comparator.hpp>
#include "antevolo_config.hpp"
#include "shm_model_utils/shm_model.hpp"
#include "mutation_strategies/no_k_neighbours.hpp"
#include "evolutionary_graph_utils/evolutionary_edge/base_evolutionary_edge.hpp"
#include "shm_model_utils/shm_model_edge_weight_calculator.hpp"

void create_console_logger() {
    using namespace logging;
    logger *lg = create_logger("");
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}

cdr_labeler::CDRLabelerConfig config;
annotation_utils::CDRAnnotatedCloneSet annotated_clone_set;

using namespace antevolo;
using namespace shm_kmer_matrix_estimator;

AntEvoloConfig antevolo_config;

class AntEvoloTest : public ::testing::Test {
public:
    void SetUp() {
        create_console_logger();
        std::string config_fname = "configs/cdr_labeler/config.info";
        cdr_labeler::load(config, config_fname);
        config.vj_finder_config.algorithm_params.germline_params.loci = "IG";
        germline_utils::GermlineDbGenerator db_generator(config.vj_finder_config.io_params.input_params.germline_input,
                                                         config.vj_finder_config.algorithm_params.germline_params);
        auto v_gene_database = db_generator.GenerateVariableDb();
        auto j_gene_database = db_generator.GenerateJoinDb();
        auto v_labeling = cdr_labeler::GermlineDbLabeler(v_gene_database, config.cdrs_params).ComputeLabeling();
        auto j_labeling = cdr_labeler::GermlineDbLabeler(j_gene_database, config.cdrs_params).ComputeLabeling();
        auto filtered_v_db = v_labeling.CreateFilteredDb();
        auto filtered_j_db = j_labeling.CreateFilteredDb();
        core::ReadArchive read_archive("test_dataset/antevolo_base_test.fastq");
        vj_finder::VJParallelProcessor processor(read_archive, config.vj_finder_config.algorithm_params,
                                                 filtered_v_db, filtered_j_db,
                                                 config.run_params.num_threads);
        auto alignment_info = processor.Process();
        config.shm_params.shm_finding_algorithm =
                cdr_labeler::CDRLabelerConfig::SHMFindingParams::SHMFindingAlgorithm::CDRFilteringSHMAlgorithm;
        cdr_labeler::ReadCDRLabeler read_labeler(config.shm_params, v_labeling, j_labeling);
        annotated_clone_set = read_labeler.CreateAnnotatedCloneSet(alignment_info);

        antevolo_config.load("configs/antevolo/config.info");
    }
};

// TEST_F(AntEvoloTest, AntEvoloBaseTest) {
//     ASSERT_EQ(annotated_clone_set.size(), 6);
//     INFO("Testing CDR3 filtered SHMs in V");
//     ASSERT_EQ(annotated_clone_set[0].VSHMs().size(), 11);
//     INFO("Testing CDR3 filtered SHMs in J");
//     ASSERT_EQ(annotated_clone_set[1].JSHMs().size(), 1);
//     INFO("Testing nested SHMs in V");
//     ASSERT_EQ(annotated_clone_set[2].VSHMs().size(), 4);
//     ASSERT_EQ(annotated_clone_set[3].VSHMs().size(), 7);
//     ASSERT_EQ(annotation_utils::SHMComparator().SHMs1AreNestedInSHMs2(annotated_clone_set[2].VSHMs(),
//                                                                       annotated_clone_set[3].VSHMs()), true);
//     INFO("Testing intersected SHMs in V");
//     ASSERT_EQ(annotation_utils::SHMComparator().GetNumberOfIntersections(annotated_clone_set[4].VSHMs(),
//                                                                          annotated_clone_set[5].VSHMs()), 15);
// }

TEST_F(AntEvoloTest, ShmModelTest) {
    INFO(antevolo_config.input_params.shm_kmer_model_igh);
    ShmModel model(antevolo_config.input_params.shm_kmer_model_igh);
    EXPECT_EQ(model.size(), 1024);

    EXPECT_EQ(model.beta_fr_params().size(), model.dirichlet_params().size());
    EXPECT_EQ(model.beta_cdr_params().size(), model.dirichlet_params().size());
    EXPECT_EQ(model.beta_full_params().size(), model.dirichlet_params().size());
    EXPECT_EQ(model.dirichlet_params().size(), model.beta_fr_success_mle().size());
    EXPECT_EQ(model.beta_fr_success_mle().size(), model.dirichlet_success_mle().size());
    EXPECT_EQ(model.beta_cdr_success_mle().size(), model.dirichlet_success_mle().size());
    EXPECT_EQ(model.beta_full_success_mle().size(), model.dirichlet_success_mle().size());
    EXPECT_EQ(model.dirichlet_success_mle().size(), model.start_point_beta_fr_params().size());
    EXPECT_EQ(model.start_point_beta_fr_params().size(), model.start_point_dirichlet_params().size());
    EXPECT_EQ(model.start_point_beta_cdr_params().size(), model.start_point_dirichlet_params().size());
    EXPECT_EQ(model.start_point_beta_full_params().size(), model.start_point_dirichlet_params().size());

    for (const auto& beta_fr_param : model.beta_fr_params()) {
        EXPECT_EQ(beta_fr_param.size(), 2);
    }

    for (const auto& beta_cdr_param : model.beta_cdr_params()) {
        EXPECT_EQ(beta_cdr_param.size(), 2);
    }

    for (const auto& beta_full_param : model.beta_full_params()) {
        EXPECT_EQ(beta_full_param.size(), 2);
    }

    for (const auto& dirichlet_param : model.dirichlet_params()) {
        EXPECT_EQ(dirichlet_param.size(), 3);
    }

    for (const auto& beta_param : model.start_point_beta_fr_params()) {
        EXPECT_EQ(beta_param.size(), 2);
    }

    for (const auto& beta_param : model.start_point_beta_cdr_params()) {
        EXPECT_EQ(beta_param.size(), 2);
    }

    for (const auto& beta_param : model.start_point_beta_full_params()) {
        EXPECT_EQ(beta_param.size(), 2);
    }

    for (const auto& dirichlet_param : model.start_point_dirichlet_params()) {
        EXPECT_EQ(dirichlet_param.size(), 3);
    }

    auto all_kmers = KmerUtils::GenerateAllKmersFixedLength(5);
    auto double_eq = [](double a, double b) {
      if (std::isnan(a)) {
          EXPECT_TRUE(std::isnan(b));
      } else {
          EXPECT_DOUBLE_EQ(a, b);
      }
    };

    for (size_t i = 0; i < model.beta_fr_params().size(); ++i) {
        double_eq(model.beta_fr_params()[i][0], model.beta_fr_params()[all_kmers[i]][0]);
        double_eq(model.beta_fr_params()[i][1], model.beta_fr_params()[all_kmers[i]][1]);

        double_eq(model.beta_cdr_params()[i][0], model.beta_cdr_params()[all_kmers[i]][0]);
        double_eq(model.beta_cdr_params()[i][1], model.beta_cdr_params()[all_kmers[i]][1]);

        double_eq(model.beta_full_params()[i][0], model.beta_full_params()[all_kmers[i]][0]);
        double_eq(model.beta_full_params()[i][1], model.beta_full_params()[all_kmers[i]][1]);
    }

    for (size_t i = 0; i < model.dirichlet_params().size(); ++i) {
        double_eq(model.dirichlet_params()[i][0], model.dirichlet_params()[all_kmers[i]][0]);
        double_eq(model.dirichlet_params()[i][1], model.dirichlet_params()[all_kmers[i]][1]);
        double_eq(model.dirichlet_params()[i][2], model.dirichlet_params()[all_kmers[i]][2]);
    }

    for (size_t i = 0; i < model.start_point_beta_fr_params().size(); ++i) {
        double_eq(model.start_point_beta_fr_params()[i][0], model.start_point_beta_fr_params()[all_kmers[i]][0]);
        double_eq(model.start_point_beta_fr_params()[i][1], model.start_point_beta_fr_params()[all_kmers[i]][1]);
        double_eq(model.start_point_beta_cdr_params()[i][0], model.start_point_beta_cdr_params()[all_kmers[i]][0]);
        double_eq(model.start_point_beta_cdr_params()[i][1], model.start_point_beta_cdr_params()[all_kmers[i]][1]);
        double_eq(model.start_point_beta_full_params()[i][0], model.start_point_beta_full_params()[all_kmers[i]][0]);
        double_eq(model.start_point_beta_full_params()[i][1], model.start_point_beta_full_params()[all_kmers[i]][1]);
    }

    for (size_t i = 0; i < model.start_point_dirichlet_params().size(); ++i) {
        double_eq(model.start_point_dirichlet_params()[i][0], model.start_point_dirichlet_params()[all_kmers[i]][0]);
        double_eq(model.start_point_dirichlet_params()[i][1], model.start_point_dirichlet_params()[all_kmers[i]][1]);
        double_eq(model.start_point_dirichlet_params()[i][2], model.start_point_dirichlet_params()[all_kmers[i]][2]);
    }
}


TEST_F(AntEvoloTest, ShmModelEdgeWeightCalculator) {
    const annotation_utils::AnnotatedClone &clone1(annotated_clone_set[2]);
    const annotation_utils::AnnotatedClone &clone2(annotated_clone_set[3]);
    // std::cout << clone1.VAlignment().query().seq << std::endl;
    // std::cout << clone2.VAlignment().query().seq << std::endl;

    BaseEvolutionaryEdge edge(clone1, clone2, 0, 0);
    ShmModel model(antevolo_config.input_params.shm_kmer_model_igh);

    AbstractMutationStrategyPtr
        mut_strategy(new NoKNeighboursMutationStrategy(antevolo_config.shm_config.mfp));
    ShmModelEdgeWeightCalculator weight_calculator(model, std::move(mut_strategy));
    double weight = weight_calculator.calculate_weigth_edge(edge);
    auto double_eq = [](double a, double b) {
      if (std::isnan(a)) {
          EXPECT_TRUE(std::isnan(b));
      } else {
          EXPECT_DOUBLE_EQ(a, b);
      }
    };
    double_eq(weight, -21.132585490945221);
}

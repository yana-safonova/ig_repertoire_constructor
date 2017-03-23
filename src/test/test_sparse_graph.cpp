#include <gtest/gtest.h>
#include <logger/log_writers.hpp>
#include "../graph_utils/sparse_graph.hpp"

struct GraphPair {
    vector<vector<bool> > matrix_;
    SparseGraphPtr sparse_;
    double fill_;

    GraphPair(size_t size, double fill) : fill_(fill) {
        for (size_t i = 0; i < size; i ++) {
            matrix_.push_back(vector<bool>(size));
        }
        for (size_t i = 0; i < size; i ++) {
            for (size_t j = 0; j < i; j ++) {
                matrix_[i][j] = matrix_[j][i] = rand() < RAND_MAX * fill_;
            }
        }
        vector<GraphEdge> edges;
        for (size_t i = 0; i < size; i ++) {
            for (size_t j = i + 1; j < size; j ++) {
                if (matrix_[i][j]) {
                    edges.push_back(GraphEdge(i, j, 1));
                }
            }
        }
        sparse_ = SparseGraphPtr(new SparseGraph(size, edges));
    }
};

class SparseGraphTestFixture : public ::testing::Test {
public:
    void SetUp();
};

vector<GraphPair> graphs;

TEST_F(SparseGraphTestFixture, TestHasEdge) {
    for (auto graph : graphs) {
        for (size_t i = 0; i < graph.matrix_.size(); i ++) {
            for (size_t j = 0; j < graph.matrix_.size(); j ++) {
                ASSERT_EQ(graph.matrix_[i][j], graph.sparse_->HasEdge(i, j));
            }
        }
    }
}

TEST_F(SparseGraphTestFixture, TestIterateEdges) {
    for (auto graph : graphs) {
        auto size = graph.matrix_.size();
        for (size_t i = 0; i < size; i ++) {
            size_t j = 0;
            for (auto v : graph.sparse_->VertexEdges(i)) {
                while (j < size && !graph.matrix_[i][j]) {
                    j ++;
                }
                ASSERT_EQ(j, v);
                j ++;
            }
            while (j < size && !graph.matrix_[i][j]) {
                j ++;
            }
            ASSERT_EQ(size, j);

            j = 0;
            for (auto itr = graph.sparse_->VertexEdges(i).begin(); itr != graph.sparse_->VertexEdges(i).end(); ++ itr) {
                while (j < size && !graph.matrix_[i][j]) {
                    j ++;
                }
                ASSERT_EQ(j, *itr);
                j ++;
            }
            while (j < size && !graph.matrix_[i][j]) {
                j ++;
            }
            ASSERT_EQ(size, j);

            j = 0;
            for (auto itr = graph.sparse_->VertexEdges(i).begin(); itr != graph.sparse_->VertexEdges(i).end(); itr ++) {
                while (j < size && !graph.matrix_[i][j]) {
                    j ++;
                }
                ASSERT_EQ(j, *itr);
                j ++;
            }
            while (j < size && !graph.matrix_[i][j]) {
                j ++;
            }
            ASSERT_EQ(size, j);
        }
    }
}

void create_console_logger() {
    using namespace logging;

    logger *lg = create_logger("");
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}

void SparseGraphTestFixture::SetUp() {
    create_console_logger();
//    time_t seed = time(NULL);
    time_t seed = 935486;
    srand(static_cast<unsigned>(seed));
    INFO("setting up with rand seed " << seed);
    vector<size_t> sizes { 0, 1, 2, 5, 10, 100 };
    for (auto size : sizes) {
        graphs.push_back(GraphPair(size, 0.0));
        double fill = 1.0;
        for (int i = 0; i < 10; i ++) {
            graphs.push_back(GraphPair(size, fill));
            fill *= 0.8;
        }
    }
}

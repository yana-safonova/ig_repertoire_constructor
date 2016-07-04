#include <verify.hpp>

#include "diversity_analyser.hpp"
#include "../graph_utils/graph_splitter.hpp"

#include <../graph_utils/graph_io.hpp>
#include <../graph_utils/graph_splitter.hpp>

namespace cdr_labeler {
    /*
    size_t DiversityAnalyser::ComputeD50(const std::vector<SparseGraphPtr> connected_components) const {
        std::vector<size_t> sizes;
        for(auto it = connected_components.begin(); it != connected_components.end(); it++)
            sizes.push_back((*it)->N());
        std::sort(sizes.begin(), sizes.end());
        size_t sum = 0;
        for(auto it = sizes.begin(); it != sizes.end(); it++) {
            sum += *it;
            if(sum >= cdr_graph_->N() / 2)
                return *it;
        }
        return size_t(-1);
    }

    void DiversityAnalyser::InitializeGraph(std::string compressed_cdr3_fasta) {
        std::string graph_fname = "cdr3.graph";
        std::string run_graph_constructor = "./build/release/bin/ig_swgraph_construct";
        std::string command_line = run_graph_constructor + " -i " + compressed_cdr3_fasta +
                " -o " + graph_fname + " --tau 3 > out";
        int err_code = system(command_line.c_str());
        VERIFY_MSG(err_code == 0, "Graph constructor finished abnormally, error code: " << err_code);
        cdr_graph_ = GraphReader(graph_fname).CreateGraph();
        ConnectedComponentGraphSplitter graph_splitter(cdr_graph_);
        auto connected_components = graph_splitter.Split();
        INFO(connected_components.size() << " connected components of CDR3 Hamming graph were created");
        size_t num_isolated = 0;
        size_t size_largest = 0;
        for(auto it = connected_components.begin(); it != connected_components.end(); it++) {
            if((*it)->N() == 1)
                num_isolated++;
            size_largest = std::max<size_t>(size_largest, (*it)->N());
        }
        INFO("# isolated CDR3s: " << num_isolated);
        INFO("Largest component size: " << size_largest);
        INFO("D50: " << ComputeD50(connected_components));
    }
     */

    void DiversityAnalyser::InitializeGraph(std::string compressed_cdr3_fasta) {
        if(compressed_cdr3_fasta == "") {
            // do something
        }
        std::string graph_fname = path::append_path(output_params_.output_dir, "cdr3_graph.graph");
        std::string run_graph_constructor = "./build/release/bin/ig_swgraph_construct";
        std::string command_line = run_graph_constructor + " -i " + compressed_cdr3_fasta +
                                   " -o " + graph_fname + " --tau 3 > " + output_params_.trash_output;
        int err_code = system(command_line.c_str());
        VERIFY_MSG(err_code == 0, "Graph constructor finished abnormally, error code: " << err_code);
        auto cdr3_graph = GraphReader(graph_fname).CreateGraph();
        ConnectedComponentGraphSplitter graph_splitter(cdr3_graph);
        cdr3_graphs_ = graph_splitter.Split();
        graph_component_map_ = cdr3_graph->GetGraphComponentMap();
        INFO(cdr3_graphs_.size() << " connected components of CDR3 Hamming graph were created");
    }

    const CompressedCDRSet& DiversityAnalyser::GetCompressedCloneSet(annotation_utils::StructuralRegion region) {
        VERIFY_MSG(region == annotation_utils::StructuralRegion::CDR1 or
                           region == annotation_utils::StructuralRegion::CDR2 or
                           region == annotation_utils::StructuralRegion::CDR3, "Region is not CDR");
        if(region == annotation_utils::StructuralRegion::CDR1)
            return cdr1_compressed_set_;
        if(region == annotation_utils::StructuralRegion::CDR2)
            return cdr2_compressed_set_;
        return cdr3_compressed_set_;
    }

    double DiversityAnalyser::ShannonIndex(annotation_utils::StructuralRegion region) {
        auto region_set = GetCompressedCloneSet(region);
        double shannon_index = 0;
        for(auto it = region_set.cbegin(); it != region_set.cend(); it++) {
            double rel_freq = double(it->second) / double(region_set.Sum());
            shannon_index += rel_freq * log(rel_freq);
        }
        return -shannon_index;
    }

    double DiversityAnalyser::SimpsonIndex(annotation_utils::StructuralRegion region) {
        auto region_set = GetCompressedCloneSet(region);
        double simpson_index = 0;
        for(auto it = region_set.cbegin(); it != region_set.cend(); it++) {
            double rel_freq = double(it->second) / double(region_set.Sum());
            simpson_index += rel_freq * rel_freq;
        }
        return simpson_index;
    }

    double DiversityAnalyser::ClonalShannonIndex() {
        double shannon_index = 0;
        for(size_t i = 0; i < cdr3_graphs_.size(); i++) {
            size_t multiplicity = 0;
            for(size_t j = 0; j < cdr3_graphs_[i]->N(); j++)
                multiplicity += cdr3_compressed_set_[graph_component_map_.GetOldVertexByNewVertex(i, j)].second;
            double rel_multiplicity = double(multiplicity) / double(cdr3_compressed_set_.Sum());
            shannon_index += rel_multiplicity * log(rel_multiplicity);
        }
        return -shannon_index;
    }

    double DiversityAnalyser::ClonalSimpsonIndex() {
        double simpson_index = 0;
        for(size_t i = 0; i < cdr3_graphs_.size(); i++) {
            size_t multiplicity = 0;
            for(size_t j = 0; j < cdr3_graphs_[i]->N(); j++)
                multiplicity += cdr3_compressed_set_[graph_component_map_.GetOldVertexByNewVertex(i, j)].second;
            double rel_multiplicity = double(multiplicity) / double(cdr3_compressed_set_.Sum());
            simpson_index += rel_multiplicity * rel_multiplicity;
        }
        return simpson_index;
    }
}
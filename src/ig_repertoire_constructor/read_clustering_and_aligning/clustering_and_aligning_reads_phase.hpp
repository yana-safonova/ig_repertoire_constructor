#pragma once

#include "stage.hpp"
#include "read_archive.hpp"
#include "ig_repertoire_constructor.hpp"
#include <fstream>

#include "aligned_read_cluster_builder.hpp"

namespace ig_repertoire_constructor {

const std::string shift_parameter_name = "shift";

class ClusteringAndAligningReadsPhase : public IgRepertoireConstructor::Phase {
public:
    ClusteringAndAligningReadsPhase() : IgRepertoireConstructor::Phase("Read clustering", "read_clustering") { }

    void run(debruijn_graph::conj_graph_pack &, const char*) {
        INFO("Clustering and aligning reads starts");

        auto read_archive_ptr = storage().GetReadArchivePtr();
        auto components_ptr = storage().GetHammingComponentsPtr();

        auto clusters_ptr = AlignedReadClusterBuilder().Build(*read_archive_ptr, *components_ptr);
        storage().SetAlignedReadClusters(clusters_ptr);

        INFO("Clustering and aligning reads ends");
    }

    void load(debruijn_graph::conj_graph_pack&,
            const std::string &load_from,
            const char* prefix) {
        std::string file_name = path::append_path(load_from, prefix) + ".aligned_clusters";
        INFO("Loading current state from " << file_name);

        auto read_archive_ptr = LoadReadArchive();
        storage().SetReadArchive(read_archive_ptr);

        std::ifstream file(file_name.c_str(), std::ios_base::in);
        VERIFY(file.good());

        VectorAlignedReadClusterPtr clusters_ptr(new std::vector <AlignedReadCluster>);
        bool need_create_new_cluster = true;
        std::string read_name;
        while (file >> read_name) {
            if (read_name == "__END__OF__CLUSTER__") {
                need_create_new_cluster = true;
            } else {
                if (need_create_new_cluster) {
                    clusters_ptr->push_back(AlignedReadCluster());
                    need_create_new_cluster = false;
                }
                std::string shiftString;
                file >> shiftString;
                int shift = atoi(shiftString.c_str() + shift_parameter_name.length() + 1); // TODO this is ugly
                clusters_ptr->back().Add(read_archive_ptr->GetReadNumberByReadName(read_name), shift);
            }
        }
        INFO(clusters_ptr->size() << " AlignedReadClusters have been loaded");
        storage().SetAlignedReadClusters(clusters_ptr);
    }

    void save(const debruijn_graph::conj_graph_pack&,
            const std::string & save_to,
            const char* prefix) const {
        std::string file_name = path::append_path(save_to, prefix) + ".aligned_clusters";
        INFO("Saving current state to " << file_name);

        auto read_archive_ptr = storage().GetReadArchivePtr();
        auto clusters_ptr = storage().GetAlignedReadClustersPtr();

        std::ofstream file;
        file.open(file_name.c_str(), std::ios_base::out);
        VERIFY(file.good());
        for (const auto & cluster : *clusters_ptr) {
            for (const auto & shifted_read : cluster) {
                file << read_archive_ptr->operator[](shifted_read.read_number).name() << "\t" << shift_parameter_name << "=" << shifted_read.shift << std::endl;
            }
            file << "__END__OF__CLUSTER__" << std::endl;
        }
        VERIFY(!file.fail());
        file.close();

        INFO(clusters_ptr->size() << " AlignedReadClusters have been saved");
    }

    virtual ~ClusteringAndAligningReadsPhase() { }
};

}

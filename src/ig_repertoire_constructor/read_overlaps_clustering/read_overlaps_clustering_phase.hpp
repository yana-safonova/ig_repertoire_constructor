#pragma once

#include "stage.hpp"
#include "read_archive.hpp"
#include "ig_repertoire_constructor.hpp"

#include "spliced_read.hpp"
#include "splicer.hpp"
#include "clusterer.hpp"

namespace ig_repertoire_constructor {

class ReadOverlapsClusteringPhase : public IgRepertoireConstructor::Phase {
public:
    ReadOverlapsClusteringPhase() :
        IgRepertoireConstructor::Phase("Build spliced reads and cluster them",
                "read_overlaps_clustering") { }

    void run(const char*) {
        ReadArchivePtr read_archive_ptr = storage().GetReadArchivePtr();

        auto clusters_ptr = storage().GetAlignedReadClustersPtr();
        std::vector <std::vector <std::vector <SplicedRead> > > spliced_reads(clusters_ptr->size());

    #pragma omp parallel for num_threads(ig_cfg::get().rp.threads_count)
        for (size_t i = 0; i < clusters_ptr->size(); ++i) {
            std::string prefix_file_name = GetPrefixFileName(i);
            Clusterer().ProcessCluster(read_archive_ptr, clusters_ptr->operator[](i), spliced_reads[i], prefix_file_name);
        }

        VectorSplicedReadClusterPtr spliced_read_clusters_ptr(new std::vector <std::vector <SplicedRead> >());
        for (const auto & components : spliced_reads) {
            spliced_read_clusters_ptr->insert(spliced_read_clusters_ptr->end(), components.begin(), components.end());
        }
        storage().SetSplicedReadClusters(spliced_read_clusters_ptr);
    }

    void load(const std::string &load_from,
            const char* prefix) {
        std::string file_name = path::append_path(load_from, prefix);
        INFO("Loading current state from " << file_name);

        auto read_archive_ptr = LoadReadArchive();
        storage().SetReadArchive(read_archive_ptr);

        auto clusters_ptr = Deserialize(file_name, read_archive_ptr);
        storage().SetSplicedReadClusters(clusters_ptr);

        INFO(clusters_ptr->size() << " clusters of SplicedReads have been loaded");
    }

    void save(const std::string & save_to,
            const char* prefix) const {
        std::string file_name = path::append_path(save_to, prefix);
        std::string file_name_numan_readable = path::append_path(save_to, prefix) + "_human_readable";
        INFO("Saving current state to " << file_name);

        auto clusters_ptr = storage().GetSplicedReadClustersPtr();
        Serialize(file_name, clusters_ptr);
        SerializeHumanReadable(file_name_numan_readable, clusters_ptr);

        INFO(clusters_ptr->size() << " clusters of SplicedReads have been saved");
    }

    virtual ~ReadOverlapsClusteringPhase() { }

private:
    static std::string GetPrefixFileName(size_t i) {
        std::ostringstream stream;
        stream << "temp_" << i << "_";
        return path::append_path(ig_cfg::get().io.temp_files, stream.str());
    }
};

}

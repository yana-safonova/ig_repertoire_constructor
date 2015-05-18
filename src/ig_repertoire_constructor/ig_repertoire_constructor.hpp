#pragma once

#include "stage.hpp"
#include "read_archive.hpp"
#include "aligned_read_cluster.hpp"
#include "spliced_read.hpp"
#include "repertoire_postprocessing/gluer.hpp"
#include "repertoire_postprocessing/aa_verification.hpp"

namespace ig_repertoire_constructor {

// place your intermediate results here
// storage is available for phases via method "storage"
class IgRepertoireConstructorStorage {
    ReadArchivePtr read_archive_ptr_;
    std::shared_ptr <std::vector <std::vector <size_t> > > hamming_components_;
    VectorAlignedReadClusterPtr vector_aligned_read_clusters_ptr_;
    VectorSplicedReadClusterPtr vector_spliced_read_clusters_ptr_;
    RepertoirePtr repertoire_ptr_;

public:
    void SetReadArchive(ReadArchivePtr read_archive_ptr) {
        read_archive_ptr_ = read_archive_ptr;
    }

    ReadArchivePtr GetReadArchivePtr() const {
        return read_archive_ptr_;
    }

    void SetHammingComponents(std::shared_ptr <std::vector <std::vector <size_t> > > hamming_components) {
        hamming_components_ = hamming_components;
    }

    std::shared_ptr <std::vector <std::vector <size_t> > > GetHammingComponentsPtr() const {
        return hamming_components_;
    }

    void SetAlignedReadClusters(VectorAlignedReadClusterPtr vector_aligned_read_clusters_ptr) {
        vector_aligned_read_clusters_ptr_ = vector_aligned_read_clusters_ptr;
    }

    VectorAlignedReadClusterPtr GetAlignedReadClustersPtr() const {
        return vector_aligned_read_clusters_ptr_;
    }

    void SetSplicedReadClusters(VectorSplicedReadClusterPtr vector_spliced_read_clusters_ptr) {
        vector_spliced_read_clusters_ptr_ = vector_spliced_read_clusters_ptr;
    }

    VectorSplicedReadClusterPtr GetSplicedReadClustersPtr() const {
        return vector_spliced_read_clusters_ptr_;
    }

    void SetRepertoirePtr(RepertoirePtr repertoire_ptr) {
        repertoire_ptr_ = repertoire_ptr;
    }

    RepertoirePtr GetRepertoirePtr() const {
        return repertoire_ptr_;
    }
};

// main stage
class IgRepertoireConstructor : public spades::CompositeStage<IgRepertoireConstructorStorage> {
public:
    IgRepertoireConstructor() : spades::CompositeStage<IgRepertoireConstructorStorage>("IgRepertoireConstructor", "ig_repertoire_constructor") { }

    void load(const std::string &,
            const char*) { }

    void save(const std::string &,
            const char*) const { }

    virtual ~IgRepertoireConstructor() { }
};

// substages



}

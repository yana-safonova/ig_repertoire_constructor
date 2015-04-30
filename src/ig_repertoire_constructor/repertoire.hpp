#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <boost/algorithm/string.hpp>

#include "sequence/sequence.hpp"
#include "io/osequencestream.hpp"
#include "read_archive.hpp"

namespace ig_repertoire_constructor {

class Cluster {
public:
    Cluster(size_t id, Sequence const &sequence) :
        id_(id),
        sequence_(sequence),
        size_(0),
        reads_(),
        reading_frame_(0) { }

    Cluster(size_t id, Sequence const &sequence, std::vector <size_t> const &reads) :
        id_(id),
        sequence_(sequence),
        size_(reads.size()),
        reads_(reads),
        reading_frame_(0) { }

    Cluster(Cluster const &cluster) :
        id_(cluster.id_),
        sequence_(cluster.sequence_),
        size_(cluster.size_),
        reads_(cluster.reads_),
        reading_frame_(cluster.reading_frame_) { }

    void AddRead(size_t read_id) {
        reads_.push_back(read_id);
        ++size_;
    }
    
    void AddReads(std::vector <size_t> const &read_ids) {
        size_ += read_ids.size();
        reads_.insert(reads_.end(), read_ids.begin(), read_ids.end());
    }

    size_t id() const {
        return id_;
    }

    size_t size() const {
        return size_;
    }

    Sequence const & sequence() const {
        return sequence_;
    }

    Sequence & sequence() {
        return sequence_;
    }

    std::vector <size_t> const &reads() const {
        return reads_;
    }

    int reading_frame() const {
        return reading_frame_;
    }

    int &reading_frame() {
        return reading_frame_;
    }

private:
    size_t id_;
    Sequence sequence_;
    size_t size_;
    std::vector <size_t> reads_;
    int reading_frame_;
};

class Repertoire {
public:
    typedef std::vector <Cluster>::iterator iterator;
    typedef std::vector <Cluster>::const_iterator const_iterator;

    Repertoire(ReadArchivePtr read_archive)  :
        read_archive_ptr_(read_archive){ }

    void AddCluster(Cluster const &cluster) {
        VERIFY_MSG(cluster_ids_to_ind_.find(cluster.id()) == cluster_ids_to_ind_.end(),
                   "Cluster with id " << cluster.id() << "already added")
        clusters_.push_back(cluster);
        cluster_ids_to_ind_[cluster.id()] = clusters_.size() - 1;
    }

    /**
     * Merge clusters with from_ind indexes to to_ind cluster.
     * Method don't delete any clusters and merged clusters from from_ind must be filtered later.
     */
    void MergeClusters(size_t to_ind, std::vector <size_t> from_ind, std::string const &new_seq, int new_reading_frame) {
        for (auto from_it = from_ind.begin(); from_it != from_ind.end(); ++from_it) {
            clusters_[to_ind].AddReads(clusters_[*from_it].reads());
        }
        clusters_[to_ind].sequence() = Sequence(new_seq.c_str());
        clusters_[to_ind].reading_frame() = new_reading_frame;
    }

    void FilterIds(std::set <size_t> const &to_filter) {
        for (auto it = clusters_.begin(); it != clusters_.end(); ) {
            if (to_filter.find(it->id()) != to_filter.end()) {
                clusters_.erase(it);
            } else {
                ++it;
            }
        }
        UpdateClusterIdsToInd();
    }

    void UpdateClusterIdsToInd() {
        cluster_ids_to_ind_.clear();
        for (size_t i = 0; i != clusters_.size(); ++i) {
            cluster_ids_to_ind_[clusters_[i].id()] = i;
        }
    }

    bool HasClusterId(size_t id) const {
        return cluster_ids_to_ind_.find(id) == cluster_ids_to_ind_.end();
    }

    Cluster const & GetClusterById(size_t id) const {
        auto cluster_id_it = cluster_ids_to_ind_.find(id);
        VERIFY_MSG(cluster_id_it != cluster_ids_to_ind_.end(), "Cluster with id " << id << " not found")
        return clusters_[cluster_id_it->second];
    }

    Cluster const & GetClusterByIndex(size_t index) const {
        return clusters_[index];
    }

    size_t size() const {
        return clusters_.size();
    }

    ReadArchivePtr read_archive_ptr() {
        return read_archive_ptr_;
    }

    iterator begin() {
        return clusters_.begin();
    }

    iterator end() {
        return clusters_.end();
    }

    const_iterator begin() const {
        return clusters_.cbegin();
    }

    const_iterator end() const {
        return clusters_.cend();
    }

    void LoadClustersFromFile(std::string const &clusters_filename);
    void LoadRcmFromFile(std::string const &rcm_filename);

    void SaveClustersToFile(std::string const &clusters_filename) const;
    void SaveRcmToFile(std::string const &rcm_filename) const;

private:
    std::vector <Cluster> clusters_;
    std::map <size_t, size_t> cluster_ids_to_ind_;
    ReadArchivePtr read_archive_ptr_;
};

typedef std::shared_ptr <Repertoire> RepertoirePtr;

void Repertoire::LoadClustersFromFile(const std::string &clusters_filename) {
    clusters_.clear();
    io::SingleStreamPtr read_stream = io::EasyStream(clusters_filename, false);
    while (!read_stream->eof()) {
        io::SingleRead read;
        (*read_stream) >> read;
        std::vector <std::string> fields;
        boost::split(fields, read.name(), boost::is_any_of("_"));
        VERIFY_MSG(fields.size() >= 10,
                   "Wrong CLUSTERS.FA file format for " << clusters_filename << ", read description must be >cluster___id___size___sz");

        size_t cluster_id = std::atoi(fields[3].c_str());
        // size_t cluster_size = std::atoi(fields[9].c_str());
        clusters_.push_back(Cluster(cluster_id, read.sequence()));
        cluster_ids_to_ind_[cluster_id] = clusters_.size() - 1;
    }
}

void Repertoire::LoadRcmFromFile(const std::string &rcm_filename) {
    std::ifstream rcm_stream(rcm_filename);
    VERIFY_MSG(rcm_stream.good(), "RCM file " << rcm_filename << " cannot be opened");
    while (!rcm_stream.eof()) {
        std::string read_name;
        size_t cluster_id;
        rcm_stream >> read_name >> cluster_id;
        if (read_name.empty()) {
            break;
        }
        size_t cluster_ind = cluster_ids_to_ind_[cluster_id];
        size_t read_id = read_archive_ptr_->GetReadNumberByReadName(read_name);
        clusters_[cluster_ind].AddRead(read_id);
    }
}

void Repertoire::SaveClustersToFile(const std::string &clusters_filename) const {
    io::osequencestream oss(clusters_filename);
    for (auto it = clusters_.begin(); it != clusters_.end(); ++it) {
        std::string cluster_name = "cluster___" + std::to_string(it->id()) + "___size___" + std::to_string(it->size());
        oss << io::SingleRead(cluster_name, it->sequence().str());
    }
}

void Repertoire::SaveRcmToFile(const std::string &rcm_filename) const {
    ofstream to_rcm(rcm_filename);
    for (auto cluster_it = clusters_.begin(); cluster_it != clusters_.end(); ++cluster_it) {
        for (auto read_it = cluster_it->reads().begin(); read_it != cluster_it->reads().end(); ++read_it) {
            to_rcm << read_archive_ptr_->operator[](*read_it).name() << '\t' << cluster_it->id() << std::endl;
        }
    }
    to_rcm.close();
}

}

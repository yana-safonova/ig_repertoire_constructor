#pragma once

#include <omp.h>
#include "sequence/sequence.hpp"
#include "spliced_read.hpp"
#include "read_archive.hpp"
#include "repertoire.hpp"
#include "upgma_clusterization.hpp"
#include "cut_vertex_clusterization.hpp"

#include "hg_clusterization.hpp"

namespace  ig_repertoire_constructor {

class PrimaryRepertoireConstructor {
    DECL_LOGGER("PrimaryRepertoireConstructor");
    VectorSplicedReadClusterPtr vector_spliced_read_clusters_ptr_;
    ReadArchivePtr read_archive_ptr_;
    size_t max_mismatches_on_tails_;
    omp_lock_t lock_;

    enum position_t {
        EMPTY,
        NO_MISMATCH,
        MISMATCH
    };

    position_t GetPositionType(int pos, std::vector <SplicedRead> const &spliced_reads) const;

    unsigned MismatchesOnLeftTail(std::vector <SplicedRead> const &spliced_reads) const;
    unsigned MismatchesOnRightTail(std::vector <SplicedRead> const &spliced_reads) const;

    bool CheckTails(std::vector <SplicedRead> const &spliced_reads) const;

    // temporary method for debug
    void PrintTracesForConsensusConstr(std::vector <int> const &occurences, std::vector <double> const &avg_qualities) const;
    char GetConsensusOnPosition(int pos, std::vector <SplicedRead> const &spliced_reads, std::vector <size_t> const &subcluster) const;
    Sequence GetConsensusSequence(std::vector <SplicedRead> const &spliced_reads, std::vector <size_t> const &subcluster) const;

public:
    PrimaryRepertoireConstructor(VectorSplicedReadClusterPtr vector_spliced_read_clusters_ptr, ReadArchivePtr read_archive_ptr, size_t max_mismatches_on_tails)
        : vector_spliced_read_clusters_ptr_(vector_spliced_read_clusters_ptr),
          read_archive_ptr_(read_archive_ptr),
          max_mismatches_on_tails_(max_mismatches_on_tails) { }

    RepertoirePtr BuildPrimaryRepertoire();
};

PrimaryRepertoireConstructor::position_t PrimaryRepertoireConstructor::GetPositionType(int pos, const std::vector<SplicedRead> &spliced_reads) const {
    bool has_reads = false;
    char prev_nucl = 0;
    for (auto it = spliced_reads.begin(); it != spliced_reads.end(); ++it) {
        if (!it->HasCharWithIndex(pos)) {
            continue;
        }
        has_reads = true;
        char curr_nucl = (*it)[pos];
        if (prev_nucl && curr_nucl != prev_nucl) {
            return MISMATCH;
        }
        prev_nucl = curr_nucl;
    }
    return (has_reads) ? NO_MISMATCH : EMPTY;
}

unsigned PrimaryRepertoireConstructor::MismatchesOnLeftTail(const std::vector<SplicedRead> &spliced_reads) const {
    unsigned mismatches = 0;
    for (int pos = -1; ; --pos) {
        position_t pos_type = GetPositionType(pos, spliced_reads);
        if (pos_type == MISMATCH) {
            ++mismatches;
        } else if (pos_type == EMPTY) {
            break;
        }
    }
    return mismatches;
}

unsigned PrimaryRepertoireConstructor::MismatchesOnRightTail(std::vector<SplicedRead> const &spliced_reads) const {
    unsigned mismatches = 0;
    for (int pos = spliced_reads[0].size(); ; ++pos) {
        position_t pos_type = GetPositionType(pos, spliced_reads);
        if (pos_type == MISMATCH) {
            ++mismatches;
        } else if (pos_type == EMPTY) {
            break;
        }
    }
    return mismatches;
}

bool PrimaryRepertoireConstructor::CheckTails(std::vector<SplicedRead> const &spliced_reads) const {
    if (MismatchesOnLeftTail(spliced_reads) + MismatchesOnRightTail(spliced_reads) <= max_mismatches_on_tails_) {
        return true;
    } else {
        return false;
    }
}

void PrimaryRepertoireConstructor::PrintTracesForConsensusConstr(std::vector <int> const &occurences,
    std::vector <double> const &avg_qualities) const {
    int zeros = 0;
    for (auto it = occurences.begin(); it != occurences.end(); ++it) {
        if (*it == 0) {
            ++zeros;
        }
    }
    if (zeros < 3) {
        TRACE("A: (" << occurences[0] << ", " << -10 * std::log10(avg_qualities[0]) << "); " <<
              "C: (" << occurences[1] << ", " << -10 * std::log10(avg_qualities[1]) << "); " <<
              "G: (" << occurences[2] << ", " << -10 * std::log10(avg_qualities[2]) << "); " <<
              "T: (" << occurences[3] << ", " << -10 * std::log10(avg_qualities[3]) << ")");
    }
}

char PrimaryRepertoireConstructor::GetConsensusOnPosition(int pos, const std::vector<SplicedRead> &spliced_reads,
                                                          std::vector <size_t> const &subcluster) const {
    std::vector <int> occurences = { 0, 0, 0, 0 };
    std::vector <double> avg_qualities = { 0, 0, 0, 0};

    size_t max_i = (subcluster.empty()) ? spliced_reads.size() : subcluster.size();
    for (size_t i = 0; i != max_i; ++i) {
        size_t spliced_read_ind = (subcluster.empty()) ? i : subcluster[i];
        SplicedRead const &curr_read = spliced_reads[spliced_read_ind];
        if (!curr_read.HasCharWithIndex(pos)) {
            continue;
        }
        char curr_nucl = curr_read[pos];
        char curr_qual = curr_read.GetQuality(pos);
        ++occurences[curr_nucl];
        avg_qualities[curr_nucl] += std::pow(10, -double(curr_qual)/10);
    }

    for (char i = 0; i != 4; ++i) {
        if (occurences[i] != 0) {
            avg_qualities[i] /= occurences[i];
        }
    }

    char best_nucl = -1;
    double best_quality = 1;
    int best_occurences = 0;
    for (char i = 0; i != 4; ++i) {
        if (occurences[i] > best_occurences) {
            best_nucl = i;
            best_quality = avg_qualities[i];
            best_occurences = occurences[i];
        } else if (occurences[i] != 0 && occurences[i] == best_occurences) {
            if (subcluster.empty()) {
                double concurrent_qual = avg_qualities[i];
                if (concurrent_qual < best_quality) {
                    best_nucl = i;
                    best_quality = avg_qualities[i];
                    best_occurences = occurences[i];
                }
            } else {
                std::vector <size_t> tmp;
                best_nucl = dignucl(GetConsensusOnPosition(pos, spliced_reads, tmp));
                best_quality = avg_qualities[i];
                best_occurences = occurences[i];
            }
        }
    }
    /*
    if (subcluster.size() > 2) {
        PrintTracesForConsensusConstr(occurences, avg_qualities);
    }
    */
    return (best_nucl == -1) ? '\0' : nucl(best_nucl);
}

Sequence PrimaryRepertoireConstructor::GetConsensusSequence(const std::vector<SplicedRead> &spliced_reads,
                                                            const std::vector<size_t> &subcluster) const {
    /*
    if (!CheckTails(spliced_reads)) {
        TRACE("Cannot construct consensus due to erroneous tails, cluster size " << spliced_reads.size());
        return Sequence("");
    }
    */
    std::string consensus_seq;

    if (spliced_reads.size() > 3) {
        std::stringstream stream;
        for (auto it = subcluster.begin(); it != subcluster.end(); ++it) {
            stream << (*it) << " ";
        }
        TRACE(stream.str());
    }

    for (int pos = GetLeftmostPosition(spliced_reads, subcluster); ; ++pos) {
        char nucl = GetConsensusOnPosition(pos, spliced_reads, subcluster);
        if (nucl == 0) {
            break;
        }
        consensus_seq.push_back(nucl);
    }

    return Sequence(consensus_seq);
}

RepertoirePtr PrimaryRepertoireConstructor::BuildPrimaryRepertoire() {
    RepertoirePtr repertoire_ptr(new Repertoire(read_archive_ptr_));
    size_t curr_cluster_id = 1;

    size_t tau = ig_cfg::get().aligning_params.overlap_mismatches_threshold;

    omp_init_lock(&lock_);
    #pragma omp parallel for num_threads(ig_cfg::get().rp.threads_count)
    for (size_t i = 0; i < vector_spliced_read_clusters_ptr_->size(); ++i) {
        std::vector <SplicedRead> const &read_group = (*vector_spliced_read_clusters_ptr_)[i];
        if (read_group.size() > 3) {
            TRACE("Process some cluster, size " << read_group.size() << ", " << read_group[0].size());
        }

        // todo: process results of cluster_constructor
        //if(i == 357134) {
        SplicedReadGroup spliced_read_group(read_group, i);
        StandardHGClustersConstructor cluster_constructor(tau, ig_cfg::get().hgc_params);
        HG_DecompositionPtr decomposition = cluster_constructor.ConstructClusters(spliced_read_group);
        //}

        continue;

        // old version
        CutVertexClusterization clusterizator(tau);
        std::shared_ptr <std::vector <std::vector <size_t> > > subclusters_ptr = clusterizator.Clusterize(read_group, i);

        for (auto subcluster_it = subclusters_ptr->begin(); subcluster_it != subclusters_ptr->end(); ++subcluster_it) {
            Sequence seq = GetConsensusSequence(read_group, *subcluster_it);
            if (seq.str().empty()) {
                TRACE("Cannot construct cluster from spliced reads of size " << subcluster_it->size());
                omp_set_lock(&lock_);
                for (auto spliced_read_ind_it = subcluster_it->begin(); spliced_read_ind_it != subcluster_it->end(); ++spliced_read_ind_it) {
                    SplicedRead const &curr_read = (read_group[*spliced_read_ind_it]);
                    repertoire_ptr->AddCluster(Cluster(curr_cluster_id,
                        (*read_archive_ptr_)[curr_read.GetReadNumber()].sequence(), {curr_read.GetReadNumber()}));
                    ++curr_cluster_id;
                }
                omp_unset_lock(&lock_);
            } else {
                std::vector <size_t> cluster_reads;
                for (auto spliced_read_ind_it = subcluster_it->begin(); spliced_read_ind_it != subcluster_it->end(); ++spliced_read_ind_it) {
                    SplicedRead const &curr_read = (read_group[*spliced_read_ind_it]);
                    cluster_reads.push_back(curr_read.GetReadNumber());
                }
                omp_set_lock(&lock_);
                Cluster cluster(curr_cluster_id, seq, cluster_reads);
                repertoire_ptr->AddCluster(cluster);
                ++curr_cluster_id;
                omp_unset_lock(&lock_);
            }
        }
    }
    assert(false);
    omp_destroy_lock(&lock_);
    return repertoire_ptr;
}


}

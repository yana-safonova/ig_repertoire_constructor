#pragma once

#include "stage.hpp"
#include "read_archive.hpp"
#include "ig_repertoire_constructor.hpp"
#include "primary_repertoire_constructor.hpp"

#define MAX_MISMATCHES_ON_TAILS 4

namespace ig_repertoire_constructor {

class PrimaryRepertoireConstructionPhase : public IgRepertoireConstructor::Phase {
public:
    PrimaryRepertoireConstructionPhase() :
        IgRepertoireConstructor::Phase("Construction of primary repertoire",
                "primary_repertoire_construction") { }

    void run(debruijn_graph::conj_graph_pack &, const char*) {
        INFO("Have " << storage().GetSplicedReadClustersPtr()->size() << " spliced reads clusters");

        PrimaryRepertoireConstructor constructor(storage().GetSplicedReadClustersPtr(),
                                                 storage().GetReadArchivePtr(), MAX_MISMATCHES_ON_TAILS);
        RepertoirePtr repertoire_ptr = constructor.BuildPrimaryRepertoire();
        storage().SetRepertoirePtr(repertoire_ptr);
        INFO(repertoire_ptr->size() << " primary clusters were constructed");
    }

    void load(debruijn_graph::conj_graph_pack&,
            const std::string &load_from,
            const char *prefix) {
        std::string prefix_file_name = path::append_path(load_from, prefix);
        INFO("Loading current_state from " << prefix_file_name);

        auto read_archive_ptr = LoadReadArchive();
        storage().SetReadArchive(read_archive_ptr);

        std::string clusters_filename = prefix_file_name + ".clusters.fa";
        std::string rcm_filename = prefix_file_name + ".rcm";
        RepertoirePtr repertoire_ptr(new Repertoire(read_archive_ptr));
        repertoire_ptr->LoadClustersFromFile(clusters_filename);
        repertoire_ptr->LoadRcmFromFile(rcm_filename);
        storage().SetRepertoirePtr(repertoire_ptr);
    }

    void save(const debruijn_graph::conj_graph_pack&,
            const std::string &save_to,
            const char *prefix) const {
        std::string prefix_file_name = path::append_path(save_to, prefix);
        INFO("Saving current state to " << prefix_file_name);
        std::string clusters_filename = prefix_file_name + ".clusters.fa";
        std::string rcm_filename = prefix_file_name + ".rcm";

        storage().GetRepertoirePtr()->SaveClustersToFile(clusters_filename);
        storage().GetRepertoirePtr()->SaveRcmToFile(rcm_filename);
    }

    virtual ~PrimaryRepertoireConstructionPhase() { }
};

}

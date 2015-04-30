#pragma once

#include "stage.hpp"
#include "read_archive.hpp"
#include "ig_repertoire_constructor.hpp"

namespace ig_repertoire_constructor {

class RepertoirePostprocessingPhase : public IgRepertoireConstructor::Phase {
public:
	RepertoirePostprocessingPhase() :
		IgRepertoireConstructor::Phase("Ig repertoire postprocessing", "repertoire_postprocessing") { }

	void run(debruijn_graph::conj_graph_pack &, const char*) {
        AAVerificator aa_verificator;
        aa_verificator.VerifyAndPredictReadingFrame(storage().GetRepertoirePtr());

        SingletonGluer::init_singleton_gluer(ig_cfg::get().sg.k, ig_cfg::get().sg.max_kmer_occurences,
                                             ig_cfg::get().sg.max_distance, ig_cfg::get().sg.min_overlap,
                                             ig_cfg::get().sg.common_mismatches_threshold, ig_cfg::get().sg.diverse_mismatches_threshold);
        SingletonGluer gluer(storage().GetRepertoirePtr(), aa_verificator.GetIncorrectClusterIds());
        gluer.GlueSingletons();

        std::string clusters_filename = path::append_path(ig_cfg::get().io.output_dir, "constructed_repertoire.clusters.fa");
        std::string rcm_filename = path::append_path(ig_cfg::get().io.output_dir, "constructed_repertoire.rcm");

        storage().GetRepertoirePtr()->SaveClustersToFile(clusters_filename);
        storage().GetRepertoirePtr()->SaveRcmToFile(rcm_filename);

        INFO("Repertoire of size " << storage().GetRepertoirePtr()->size() << " was constructed");
        INFO("CLUSTERS.FA file was written to " << clusters_filename);
        INFO("RCM file was written to " << rcm_filename);
    }

	void load(debruijn_graph::conj_graph_pack&,
            const std::string &/*load_from*/,
            const char */*prefix*/) {
        /*
        std::string file_name = path::append_path(load_from, prefix);
        INFO("Loading current state from " << file_name);
        std::string clusters_filename = file_name + ".clusters.fa";
        std::string rcm_filename = file_name + ".rcm";

        auto read_archive_ptr = LoadReadArchive();

        RepertoirePtr repertoire_ptr(new Repertoire(read_archive_ptr));
        repertoire_ptr->LoadClustersFromFile(clusters_filename);
        repertoire_ptr->LoadRcmFromFile(rcm_filename);
        storage().SetRepertoirePtr(repertoire_ptr);
        */
    }

	void save(const debruijn_graph::conj_graph_pack&,
            const std::string &/*save_to*/,
            const char* /*prefix*/) const {
        /*
        std::string file_name = path::append_path(save_to, prefix);
        INFO("Saving current state to " << file_name);
        std::string clusters_filename = file_name + ".clusters.fa";
        std::string rcm_filename = file_name + ".rcm";

        storage().GetRepertoirePtr()->SaveClustersToFile(clusters_filename);
        storage().GetRepertoirePtr()->SaveRcmToFile(rcm_filename);
        */
    }

	virtual ~RepertoirePostprocessingPhase() { }
};

}


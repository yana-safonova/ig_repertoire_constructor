#pragma once

#include "ig_repertoire_constructor.hpp"
#include "hamming_graph_building/hamming_graph_building_phase.hpp"
#include "read_clustering_and_aligning/clustering_and_aligning_reads_phase.hpp"
#include "read_overlaps_clustering/read_overlaps_clustering_phase.hpp"
#include "primary_repertoire_construction/primary_repertoire_construction_phase.hpp"
#include "repertoire_postprocessing/repertoire_postprocessing_phase.hpp"

namespace ig_repertoire_constructor {

void run_ig_repertoire_constructor() {
    INFO("IgRepertoireConstructor started");
    spades::StageManager Ig_Manager ( {
    	ig_cfg::get().rp.developer_mode,
    	ig_cfg::get().io.load_from,
    	ig_cfg::get().io.output_saves } );
    auto main_ig_phase = new IgRepertoireConstructor();
    main_ig_phase->
        add(new HammingGraphBuildingPhase()) ->
    	add(new ClusteringAndAligningReadsPhase()) ->
    	add(new ReadOverlapsClusteringPhase()) ->
    	add(new PrimaryRepertoireConstructionPhase())
    	//-> add(new RepertoirePostprocessingPhase())
    	;

    Ig_Manager.add(main_ig_phase);

    Ig_Manager.run(ig_cfg::get().rp.entry_point.c_str());
    INFO("IgRepertoireConstructor finished");
}

}

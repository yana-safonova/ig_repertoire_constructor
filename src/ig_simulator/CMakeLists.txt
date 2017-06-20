project(ig_simulator CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${SPADES_MAIN_INCLUDE_DIR})
include_directories(${CORE_DIR})
include_directories(${VDJ_UTILS_DIR})
include_directories(${VJ_FINDER_DIR})
include_directories(${CDR_LABELER_DIR})

add_library(ig_simulator_library STATIC
        ig_simulator_config.cpp
        ig_simulator_launch.cpp
        base_repertoire/gene_chooser/abstract_gene_chooser.cpp
        base_repertoire/gene_chooser/uniform_gene_chooser.cpp
        base_repertoire/gene_chooser/config_based_getter.cpp
        base_repertoire/nucleotides_remover/abstract_nucleotides_remover.cpp
        base_repertoire/nucleotides_remover/uniform_nucleotides_remover.cpp
        base_repertoire/nucleotides_remover/config_based_getter.cpp
        simulation_routines.cpp
        base_repertoire/p_nucleotides_creator/abstract_nucleotides_creator.cpp
        base_repertoire/p_nucleotides_creator/uniform_nucleotides_creator.cpp
        base_repertoire/p_nucleotides_creator/config_based_getter.cpp
        base_repertoire/metaroot_creator/metaroot_creator.cpp
        base_repertoire/metaroot_creator/config_based_getter.cpp
        base_repertoire/n_nucleotides_inserter/abstract_n_nucleotides_inserter.cpp
        base_repertoire/n_nucleotides_inserter/uniform_n_nucleotides_inserter.cpp
        base_repertoire/n_nucleotides_inserter/config_based_getter.cpp
        base_repertoire/metaroot/metaroot.cpp
        base_repertoire/productivity_checker/productivity_checker.cpp
        base_repertoire/multiplicity_creator/multiplicity_creator.cpp
        base_repertoire/base_repertoire_simulator.cpp
        base_repertoire/metaroot_cluster/metaroot_cluster.cpp
        base_repertoire/base_repertoire.cpp
        clonal_trees/tree/node.cpp
        clonal_trees/tree_creator/pool_manager.cpp
        clonal_trees/tree_creator/cartesian_tree.cpp
        clonal_trees/tree/tree.cpp
        clonal_trees/forest/forest.cpp
        clonal_trees/tree_creator/tree_creator.cpp
        clonal_trees/tree_creator/tree_size_generator.cpp
        clonal_trees/tree_creator/shm_creator.cpp
        clonal_trees/tree_creator/forest_creator.cpp
        clonal_trees/tree_creator/forest_storage_creator.cpp
        clonal_trees/tree_creator/exporters.cpp
        clonal_trees/fast_stop_codon_checker/fast_stop_codon_checker.cpp
        base_repertoire/gene_chooser/custom_gene_chooser.cpp)

target_link_libraries(ig_simulator_library
        vj_finder_library
        cdr_labeler_library
        input
        ${COMMON_LIBRARIES}
        )


add_executable(ig_simulator main.cpp)

target_link_libraries(ig_simulator ig_simulator_library)

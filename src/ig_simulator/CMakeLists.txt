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
        gene_chooser/abstract_gene_chooser.cpp
        gene_chooser/uniform_gene_chooser.cpp
        nucleotides_remover/abstract_nucleotides_remover.cpp
        nucleotides_remover/uniform_nucleotides_remover.cpp
        random_index.cpp
        p_nucleotides_creator/abstract_nucleotides_creator.cpp
        p_nucleotides_creator/uniform_nucleotides_creator.cpp
        metaroot_creator/metaroot_creator.cpp
        n_nucleotides_inserter/abstract_n_nucleotides_inserter.cpp
        n_nucleotides_inserter/uniform_n_nucleotides_inserter.cpp
        metaroot/metaroot.cpp)

target_link_libraries(ig_simulator_library
        vj_finder_library
        cdr_labeler_library
        input
        ${COMMON_LIBRARIES}
        )


add_executable(ig_simulator main.cpp)

target_link_libraries(ig_simulator ig_simulator_library)

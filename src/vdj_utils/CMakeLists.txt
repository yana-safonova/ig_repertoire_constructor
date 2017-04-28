project(vdj_utils CXX)

INCLUDE_DIRECTORIES(SYSTEM ${SEQAN_INCLUDE_DIR})
include_directories(${IGREC_MAIN_INCLUDE_DIR})
include_directories(${CORE_DIR})

add_library(vdj_utils STATIC
        germline_utils/lymphocyte_type.cpp
        germline_utils/chain_type.cpp
        germline_utils/germline_gene_type.cpp
        germline_utils/germline_db_generator.cpp
        germline_utils/germline_databases/immune_gene_database.cpp
        germline_utils/germline_databases/chain_database.cpp
        germline_utils/germline_databases/custom_gene_database.cpp
        alignment_utils/alignment_positions.cpp
        alignment_utils/pairwise_alignment.cpp
        annotation_utils/cdr_labeling_primitives.cpp
        annotation_utils/shm_annotation/shm_annotation.cpp
        annotation_utils/shm_annotation/shm_calculator.cpp
        annotation_utils/aa_annotation/aa_annotation.cpp
        annotation_utils/aa_annotation/aa_calculator.cpp
        annotation_utils/annotated_clone.cpp
        annotation_utils/annotated_clone_calculator.cpp
        annotation_utils/annotated_clone_set.cpp
        germline_utils/germline_config.cpp)

target_link_libraries(vdj_utils core ${COMMON_LIBRARIES})

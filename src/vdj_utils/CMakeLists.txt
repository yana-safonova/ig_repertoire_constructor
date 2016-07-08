project(vdj_utils CXX)

INCLUDE_DIRECTORIES(SYSTEM ${SEQAN_INCLUDE_DIR})
include_directories(${SPADES_MAIN_INCLUDE_DIR})

add_library(vdj_utils STATIC
        germline_utils/lymphocyte_type.cpp
        germline_utils/chain_type.cpp
        germline_utils/germline_gene_type.cpp
        germline_utils/germline_databases/immune_gene_database.cpp
        germline_utils/germline_databases/chain_database.cpp
        germline_utils/germline_databases/custom_gene_database.cpp
        )
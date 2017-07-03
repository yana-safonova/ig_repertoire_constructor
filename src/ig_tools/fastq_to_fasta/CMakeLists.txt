project(fastq_to_fasta CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${IGREC_MAIN_INCLUDE_DIR})
include_directories(${IGREC_MAIN_SRC_DIR}/ig_tools)

add_executable(fastq_to_fasta
               fastq_to_fasta.cpp)

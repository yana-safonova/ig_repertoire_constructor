project(umi_experiments CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_executable(umi_to_fastq umi_to_fastq.cpp)
project(umi_experiments CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

link_libraries(
        input
        ${COMMON_LIBRARIES}
        graph_utils
        boost_filesystem
        boost_system
        boost_program_options
)

# needed for CLion to include them into the project
file(GLOB HEADER_FILES **/*.hpp)

add_executable(analyze_umi_groups analyze_umi_groups.cpp ../umi_experiments/utils.cpp ../umi_experiments/umi_utils.cpp ${HEADER_FILES} ion_utils.hpp)
add_executable(compress_ion compress_ion.cpp ../umi_experiments/utils.cpp ../umi_experiments/umi_utils.cpp ${HEADER_FILES} ion_utils.hpp)
add_executable(compress_germline tools/compress_germline.cpp ../umi_experiments/utils.cpp ion_utils.hpp)
add_executable(uncompress_vjf_results uncompress_vjf_results.cpp ../umi_experiments/utils.cpp ../umi_experiments/utils/io.cpp)

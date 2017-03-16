project(ig_simulator CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${SPADES_MAIN_INCLUDE_DIR})
include_directories(${CORE_DIR})
include_directories(${VDJ_UTILS_DIR})
include_directories(${VJ_FINDER_DIR})

add_executable(ig_simulator
        ig_simulator_config.cpp
        ig_simulator_launch.cpp
        main.cpp)

target_link_libraries(ig_simulator
        vj_finder_library
        input
        ${COMMON_LIBRARIES}
        )

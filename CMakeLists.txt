project(ig_repertoire_constructor)
cmake_minimum_required(VERSION 2.8)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RelWithAsserts" CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE)

add_subdirectory(configs)
add_subdirectory(src)

if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "/opt" CACHE PATH "default install path" FORCE)
endif()
set(INSTALL_DIR "igrec")

install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/VERSION"
              "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE"
              "${CMAKE_CURRENT_SOURCE_DIR}/README"
              "${CMAKE_CURRENT_SOURCE_DIR}/GPLv2.txt"
              "${CMAKE_CURRENT_SOURCE_DIR}/init.py"
              "${CMAKE_CURRENT_SOURCE_DIR}/manual.html"
        DESTINATION "${INSTALL_DIR}"
        COMPONENT runtime)

install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/test_dataset"
                  "${CMAKE_CURRENT_SOURCE_DIR}/data"
                  "${CMAKE_CURRENT_SOURCE_DIR}/pipeline_makefiles"
                  "${CMAKE_CURRENT_SOURCE_DIR}/docs"
                  "${CMAKE_CURRENT_SOURCE_DIR}/configs"
        DESTINATION "${INSTALL_DIR}"
        COMPONENT runtime)

install(DIRECTORY "${CMAKE_BINARY_DIR}/bin"
        DESTINATION "${INSTALL_DIR}/build/release"
        COMPONENT runtime
		USE_SOURCE_PERMISSIONS)

install(PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/igrec.py"
                 "${CMAKE_CURRENT_SOURCE_DIR}/mass_spectra_analyzer.py"
                 "${CMAKE_CURRENT_SOURCE_DIR}/dense_subgraph_finder.py"
                 "${CMAKE_CURRENT_SOURCE_DIR}/igquast.py"
                 "${CMAKE_CURRENT_SOURCE_DIR}/barcoded_igrec.py"
                 "${CMAKE_CURRENT_SOURCE_DIR}/diversity_analyzer.py"
        DESTINATION "${INSTALL_DIR}"
        COMPONENT runtime)

install(PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/src/umi_experiments/stats/print_umi_graph_stats.py"
        DESTINATION "${INSTALL_DIR}/src/umi_experiments/stats/"
        COMPONENT runtime)

install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/py"
        DESTINATION "${INSTALL_DIR}"
        COMPONENT runtime
		USE_SOURCE_PERMISSIONS)

# Make symlink if possible (return 0 in any case)
install(CODE "EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E create_symlink \"\${CMAKE_INSTALL_PREFIX}/${INSTALL_DIR}/igrec.py\" \"/usr/bin/igrec.py\")")

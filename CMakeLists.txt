project(igrec)
cmake_minimum_required(VERSION 2.8)

include(ExternalProject)

ExternalProject_Add(gperftools-package
                    gperftools
                    URL https://github.com/gperftools/gperftools/releases/download/gperftools-2.6.1/gperftools-2.6.1.tar.gz
                    DOWNLOAD_DIR ${CMAKE_BINARY_DIR}/thirdparty-dist
                    INSTALL_DIR ${CMAKE_BINARY_DIR}/thirdparty
                    URL_MD5 e010d656b8b62e61c0af0c4d79843782
                    CONFIGURE_COMMAND ./configure --enable-shared=yes --enable-static=no --prefix=<INSTALL_DIR>
                    BUILD_COMMAND ${CMAKE_MAKE_PROGRAM}
                    BUILD_IN_SOURCE 1
                    INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install && ${CMAKE_MAKE_PROGRAM} clean
                    LOG_DOWNLOAD 1
                    LOG_BUILD 1)

add_library(profiler SHARED IMPORTED GLOBAL)
# If you want to automatically enable profiling, you need link against SHARED library. See https://github.com/gperftools/gperftools/issues/415
# GLOBAL means available in ALL CMakeLists included in this one
set_property(TARGET profiler PROPERTY IMPORTED_LOCATION "${CMAKE_BINARY_DIR}/thirdparty/lib/${CMAKE_SHARED_LIBRARY_PREFIX}profiler${CMAKE_SHARED_LIBRARY_SUFFIX}")
add_dependencies(profiler gperftools-package)

# set_target_properties(profiler PROPERTIES EXCLUDE_FROM_ALL 1 EXCLUDE_FROM_DEFAULT_BUILD 1)
set_target_properties(gperftools-package PROPERTIES EXCLUDE_FROM_ALL 1 EXCLUDE_FROM_DEFAULT_BUILD 1)

if(CMAKE_GOOGLE_PROFILER)
  link_libraries(profiler)  # It is necessary for correct dependency working
  get_target_property(PROFILER profiler LOCATION)
  link_libraries("-Wl,--no-as-needed,${PROFILER},--as-needed")  # Force linking (otherwise linker may exclude profiler lib automatically and profiling will not work)
endif()


ExternalProject_Add(valgrind-package
                    PREFIX valgrind
                    URL http://valgrind.org/downloads/valgrind-3.12.0.tar.bz2
                    DOWNLOAD_DIR ${CMAKE_BINARY_DIR}/thirdparty-dist
                    INSTALL_DIR ${CMAKE_BINARY_DIR}/thirdparty
                    URL_MD5 6eb03c0c10ea917013a7622e483d61bb
                    CONFIGURE_COMMAND ./configure --prefix=<INSTALL_DIR>
                    BUILD_COMMAND ${CMAKE_MAKE_PROGRAM}
                    BUILD_IN_SOURCE 1
                    INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install && ${CMAKE_MAKE_PROGRAM} clean
                    LOG_DOWNLOAD 1
                    LOG_BUILD 1)

add_executable(valgrind IMPORTED GLOBAL)
set_property(TARGET valgrind PROPERTY IMPORTED_LOCATION "${CMAKE_BINARY_DIR}/thirdparty/bin/valgrind${CMAKE_EXECUTABLE_SUFFIX}")
add_dependencies(valgrind valgrind-package)
set_target_properties(valgrind PROPERTIES EXCLUDE_FROM_ALL 1 EXCLUDE_FROM_DEFAULT_BUILD 1)
set_target_properties(valgrind-package PROPERTIES EXCLUDE_FROM_ALL 1 EXCLUDE_FROM_DEFAULT_BUILD 1)


if(CMAKE_NVERIFY)
  MESSAGE("VERIFYs disabled")
  add_definitions(-DNVERIFY)
endif(CMAKE_NVERIFY)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RelWithAsserts" CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE)

add_subdirectory(configs)
add_subdirectory(src)

if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "/usr/local" CACHE PATH "default install path" FORCE)
endif()

# Expand tildas and convert to absolute path
FILE(TO_CMAKE_PATH "${CMAKE_INSTALL_PREFIX}" CMAKE_INSTALL_PREFIX)

set(INSTALL_DIR "share/igrec")

install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/VERSION"
              "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE"
              "${CMAKE_CURRENT_SOURCE_DIR}/README"
              "${CMAKE_CURRENT_SOURCE_DIR}/GPLv2.txt"
              "${CMAKE_CURRENT_SOURCE_DIR}/init.py"
        DESTINATION "${INSTALL_DIR}"
        COMPONENT runtime)

install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/manual.html"
              "${CMAKE_CURRENT_SOURCE_DIR}/diversity_analyzer_manual.html"
              "${CMAKE_CURRENT_SOURCE_DIR}/dsf_manual.html"
              "${CMAKE_CURRENT_SOURCE_DIR}/barcodedIgReC_manual.html"
              "${CMAKE_CURRENT_SOURCE_DIR}/igquast_manual.html"
        DESTINATION "${INSTALL_DIR}"
        COMPONENT manual)

install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/docs"
        DESTINATION "${INSTALL_DIR}"
        COMPONENT manual
        PATTERN "*.pdf" EXCLUDE)

install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/test_dataset"
                  "${CMAKE_CURRENT_SOURCE_DIR}/data"
                  "${CMAKE_CURRENT_SOURCE_DIR}/pipeline_makefiles"
                  "${CMAKE_CURRENT_SOURCE_DIR}/configs"
        DESTINATION "${INSTALL_DIR}"
        COMPONENT runtime
        PATTERN "*.template" EXCLUDE
        PATTERN "CMakeLists.txt" EXCLUDE
        PATTERN ".gitignore" EXCLUDE)

install(DIRECTORY "${CMAKE_BINARY_DIR}/bin"
        DESTINATION "${INSTALL_DIR}/build/release"
        COMPONENT runtime
        USE_SOURCE_PERMISSIONS
        PATTERN "test_*" EXCLUDE)

install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/py"
        DESTINATION "${INSTALL_DIR}"
        COMPONENT runtime
        USE_SOURCE_PERMISSIONS
        PATTERN "*.pyc" EXCLUDE
        PATTERN ".gitignore" EXCLUDE)

set(tools "igrec.py;mass_spectra_analyzer.py;dense_subgraph_finder.py;igquast.py;barcoded_igrec.py;diversity_analyzer.py;ig_simulator.py;antevolo.py;tree_visualizer.py;igrec_gui.py")
foreach(tool ${tools})
    install(PROGRAMS ${tool}
            DESTINATION "${INSTALL_DIR}"
            COMPONENT runtime)
endforeach()

foreach(tool ${tools})
    # Conctruct fake symlinks
    EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E create_symlink "../share/igrec/${tool}" "${CMAKE_BINARY_DIR}/${tool}")
    # ...and add them into the installation
    # I know, I know, I definitely know how much dirty it is
    # Unfortunately, cpack cannot trace symlinks automatically and cmake-install does not support smth like install(SYMLINK)
    install(FILES "${CMAKE_BINARY_DIR}/${tool}"
            DESTINATION "bin"
            COMPONENT runtime)
endforeach()

# Disabled. We have already created symlinks manually above
# foreach(tool ${tools})
#     install(CODE "EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E create_symlink \"${CMAKE_INSTALL_PREFIX}/${INSTALL_DIR}/${tool}\" \"${CMAKE_INSTALL_PREFIX}/bin/${tool}\")")
# endforeach()

set(CPACK_PACKAGE_NAME "igrec")
set(CPACK_PACKAGING_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")
set(CPACK_GENERATOR "TGZ;ZIP;DEB;RPM")
set(CPACK_SOURCE_GENERATOR "")
set(CPACK_OUTPUT_FILE_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/packages")  # TODO add an arg
set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/README")

FILE(READ "${CMAKE_CURRENT_SOURCE_DIR}/VERSION" VERSION)
STRING(STRIP "${VERSION}" VERSION)
set(CPACK_PACKAGE_VERSION ${VERSION})
# Are they really necessary?
# SET(CPACK_PACKAGE_VERSION_MAJOR "3")
# SET(CPACK_PACKAGE_VERSION_MINOR "0")
# SET(CPACK_PACKAGE_VERSION_PATCH "0")

set(CPACK_DEBIAN_PACKAGE_DEPENDS "python")
# TODO Add seaborn somehow
# set(CPACK_DEBIAN_PACKAGE_SUGGESTS "python-matplotlib, python-biopython, python-numpy, python-scipy, python-pandas")

set(CPACK_PACKAGE_CONTACT "Alexander Shlemov <shlemovalex@gmail.com>")
set(CPACK_DEBIAN_PACKAGE_HOMEPAGE "http://yana-safonova.github.io/ig_repertoire_constructor/")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Y-tools: tools for construction of full-length adaptive immune repertoires and immunoinformatics analysis")
# TODO: Add full description
set(CPACK_PACKAGE_DESCRIPTION "Y-tools: tools for construction of full-length adaptive immune repertoires and immunoinformatics analysis")

set(CPACK_PACKAGE_VENDOR "Center for Algorithmic Biotechnology, Saint Petersburg University, Russia")

# trace debug information during CPackDeb run
set(CPACK_DEBIAN_PACKAGE_DEBUG TRUE) 

# RPM
# trace debug information during CPackRPM run
set(CPACK_RPM_PACKAGE_LICENSE "GPLv2")
set(CPACK_RPM_PACKAGE_DEBUG TRUE)
set(CPACK_RPM_PACKAGE_DESCRIPTION ${CPACK_PACKAGE_DESCRIPTION})
set(CPACK_RPM_PACKAGE_REQUIRES ${CPACK_DEBIAN_PACKAGE_DEPENDS}) # TODO check
# set(CPACK_RPM_PACKAGE_SUGGESTS ${CPACK_DEBIAN_PACKAGE_SUGGESTS})
set(CPACK_RPM_PACKAGE_URL ${CPACK_DEBIAN_PACKAGE_HOMEPAGE})


set(CPACK_STRIP_FILES TRUE)
set(CPACK_PACKAGE_EXECUTABLES $tools)

set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
set(CPACK_MONOLITHIC_INSTALL TRUE)

# Package parameters should be set BEFORE include(CPack)
include(CPack)

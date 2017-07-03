# -*- cmake -*-

# Default configuration
set(IGREC_DEFAULT_BUILD_TYPE "RelWithAsserts" CACHE STRING "IgReC default build type")
if (NOT CMAKE_BUILD_TYPE)
  message("Setting default build configuration: ${IGREC_DEFAULT_BUILD_TYPE}")
  set(CMAKE_BUILD_TYPE "${IGREC_DEFAULT_BUILD_TYPE}" CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithAsserts RelWithDebInfo."
      FORCE)
endif()

# Define option for turning on/off debug logging
option(IGREC_DEBUG_LOGGING "Turn on debug / trace logging" ON)

# Define option for static / dynamic build.
option(IGREC_STATIC_BUILD "Link IgReC statically" OFF)
if (IGREC_STATIC_BUILD)
  # it'll make cmake to find libraries archives, not dynamic link
  set(CMAKE_FIND_LIBRARY_SUFFIXES .a) 
  set(LINK_SEARCH_START_STATIC TRUE)
  set(LINK_SEARCH_END_STATIC TRUE)

  if (APPLE)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc")
  else()
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
    add_definitions(-static)
  endif()

  set(Boost_USE_STATIC_LIBS        ON)
  set(Boost_USE_STATIC_RUNTIME     ON)
endif()

# Define minimum and maximum K
set(IGREC_MIN_K 1 CACHE INTEGER "Minimum k-mer length")
set(IGREC_MAX_K 128 CACHE INTEGER "Maximum k-mer length")
configure_file("${IGREC_MAIN_INCLUDE_DIR}/k_range.hpp.in"
               "${IGREC_BUILT_INCLUDE_DIR}/k_range.hpp")

# Define boost root
set(IGREC_BOOST_ROOT "" CACHE PATH "Boost root used to build IgReC")

# Various internal stuff
option(IGREC_BUILD_INTERNAL "Build internal projects" OFF)
option(IGREC_USE_TCMALLOC "Link IgReC with TCMalloc" OFF)
if (IGREC_USE_TCMALLOC)
  find_package(GooglePerfTools REQUIRED)

  if (GOOGLE_PERFTOOLS_ROOT)
    # add the automatically determined parts of the RPATH
    # which point to directories outside the build tree to the install RPATH
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  endif()
endif()
option(IGREC_USE_JEMALLOC "Link IgReC with JEMalloc" ON)
if (IGREC_USE_JEMALLOC AND IGREC_USE_TCMALLOC)
  message(FATAL_ERROR "Cannot link to both TCMalloc and JEMalloc")
endif()

############################################################################
# Copyright (c) 2015-2015 Saint-Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(vdj_labeler_0.3 CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${SPADES_MAIN_INCLUDE_DIR})

find_path(SEQAN_INCLUDE_DIR seqan/seq_io.h
  HINTS ../../ext/include ${EXT_DIR}/include)

if(SEQAN_INCLUDE_DIR)
  message("seqan headers found in ${SEQAN_INCLUDE_DIR}")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -Wextra")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

INCLUDE_DIRECTORIES(SYSTEM ${SEQAN_INCLUDE_DIR})

add_executable(vdj_labeler_0.3
               vdj_config.cpp
               ../vj_finder/vj_finder_config.cpp
               main.cpp)

target_link_libraries(vdj_labeler_0.3 input yaml-cpp ${COMMON_LIBRARIES})

if (SPADES_STATIC_BUILD)
  set_target_properties(vdj_labeler_0.3 PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

install(TARGETS vdj_labeler_0.3
        DESTINATION bin
        COMPONENT runtime)
install(DIRECTORY "${SPADES_CFG_DIR}/vdj_labeler_0.3"
        DESTINATION share/spades/configs
        FILES_MATCHING PATTERN "*.info.template")

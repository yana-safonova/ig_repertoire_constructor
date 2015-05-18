############################################################################
# Copyright (c) 2011-2014 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(ig_repertoire_constructor CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${SPADES_MAIN_INCLUDE_DIR})

add_executable(ig_repertoire_constructor
               ${SPADES_MAIN_INCLUDE_DIR}/stage.cpp
               read_archive.cpp
               ig_config.cpp          
               main.cpp)

target_link_libraries(ig_repertoire_constructor input mph_index BamTools yaml-cpp nlopt ${COMMON_LIBRARIES} ssw)

if (SPADES_STATIC_BUILD)
  set_target_properties(ig_repertoire_constructor PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

install(TARGETS ig_repertoire_constructor
        DESTINATION bin
        COMPONENT runtime)
install(DIRECTORY "${SPADES_CFG_DIR}/ig_repertoire_constructor"
        DESTINATION share/spades/configs
        FILES_MATCHING PATTERN "*.info.template")

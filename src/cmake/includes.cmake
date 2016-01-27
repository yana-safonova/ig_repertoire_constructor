# -*- cmake -*-

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_INCLUDE_SYSTEM_FLAG_C "-isystem ")
set(CMAKE_INCLUDE_SYSTEM_FLAG_CXX "-isystem ")
include_directories(${SPADES_MAIN_INCLUDE_DIR} ${SPADES_BUILT_INCLUDE_DIR})
include_directories(SYSTEM "${EXT_DIR}/include")

if (SPADES_USE_TCMALLOC)
  include_directories("${GOOGLE_PERFTOOLS_INCLUDE_DIR}")
endif()

project(build_info)
cmake_minimum_required(VERSION 2.8)


include(rpavlik-cmake-modules/GetGitRevisionDescription.cmake)
git_describe(BUILD_INFO_GIT_DESCRIBE)
get_git_head_revision(BUILD_INFO_GIT_REFSPEC BUILD_INFO_GIT_HASH)
git_get_exact_tag(BUILD_INFO_GIT_TAG)
STRING(SUBSTRING ${BUILD_INFO_GIT_HASH} 0 7 BUILD_INFO_GIT_HASH7)

SET(BUILD_INFO_CMAKE ${CMAKE_VERSION})
SET(BUILD_INFO_SYSTEM ${CMAKE_SYSTEM})
SET(BUILD_INFO_CPU ${CMAKE_SYSTEM_PROCESSOR})
SET(BUILD_INFO_C_COMPILER_ID ${CMAKE_C_COMPILER_ID})
SET(BUILD_INFO_CXX_COMPILER_ID ${CMAKE_CXX_COMPILER_ID})
SET(BUILD_INFO_C_COMPILER ${CMAKE_C_COMPILER})
SET(BUILD_INFO_CXX_COMPILER ${CMAKE_CXX_COMPILER})
SET(BUILD_INFO_BUILD_TYPE ${CMAKE_BUILD_TYPE})
SET(BUILD_INFO_C_FLAGS ${CMAKE_C_FLAGS})
SET(BUILD_INFO_CXX_FLAGS ${CMAKE_CXX_FLAGS})

if(EXISTS "${CMAKE_SOURCE_DIR}/../VERSION")
  FILE(READ "${CMAKE_SOURCE_DIR}/../VERSION" BUILD_INFO_VERSION)
  STRING(STRIP ${BUILD_INFO_VERSION} BUILD_INFO_VERSION)
else()
  SET(BUILD_INFO_VERSION "<unknown>")
endif()


configure_file("${CMAKE_CURRENT_SOURCE_DIR}/build_info.py.in"
               "${IGREC_MAIN_SRC_DIR}/../py/build_info/build_info.py"
               @ONLY)


configure_file("${CMAKE_CURRENT_SOURCE_DIR}/build_info.cpp.in"
               "${CMAKE_CURRENT_SOURCE_DIR}/build_info.cpp"
               @ONLY)
add_library(build_info STATIC
            build_info.cpp)


add_executable(test_build_info test_build_info.cpp)
target_link_libraries(test_build_info build_info)

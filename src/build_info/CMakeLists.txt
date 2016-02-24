cmake_minimum_required(VERSION 2.8)

include(rpavlik-cmake-modules/GetGitRevisionDescription.cmake)


git_describe(BUILD_INFO_GIT_DESCRIBE)
get_git_head_revision(BUILD_INFO_GIT_REFSPEC BUILD_INFO_GIT_HASH)
git_get_exact_tag(BUILD_INFO_GIT_TAG)
STRING(SUBSTRING ${BUILD_INFO_GIT_HASH} 0 7 BUILD_INFO_GIT_HASH7)

# message("${BUILD_INFO_GIT_HASH}")
# message("${BUILD_INFO_GIT_REFSPEC}")
# message("${BUILD_INFO_GIT_DESCRIBE}")
# message("${BUILD_INFO_GIT_TAG}")

add_definitions(-D BUILD_INFO_GIT_HASH='"${BUILD_INFO_GIT_HASH}"')
add_definitions(-D BUILD_INFO_GIT_HASH7='"${BUILD_INFO_GIT_HASH7}"')
add_definitions(-D BUILD_INFO_GIT_REFSPEC='"${BUILD_INFO_GIT_REFSPEC}"')
add_definitions(-D BUILD_INFO_GIT_DESCRIBE='"${BUILD_INFO_GIT_DESCRIBE}"')
add_definitions(-D BUILD_INFO_GIT_TAG='"${BUILD_INFO_GIT_TAG}"')

add_definitions(-D BUILD_INFO_CMAKE='"${CMAKE_VERSION}"')
add_definitions(-D BUILD_INFO_SYSTEM='"${CMAKE_SYSTEM}"')
add_definitions(-D BUILD_INFO_CPU='"${CMAKE_SYSTEM_PROCESSOR}"')

add_definitions(-D BUILD_INFO_C_COMPILER_ID='"${CMAKE_C_COMPILER_ID}"')
add_definitions(-D BUILD_INFO_CXX_COMPILER_ID='"${CMAKE_CXX_COMPILER_ID}"')
add_definitions(-D BUILD_INFO_C_COMPILER='"${CMAKE_C_COMPILER}"')
add_definitions(-D BUILD_INFO_CXX_COMPILER='"${CMAKE_CXX_COMPILER}"')

add_definitions(-D BUILD_INFO_BUILD_TYPE='"${CMAKE_BUILD_TYPE}"')
SET(BUILD_INFO_C_FLAGS ${CMAKE_C_FLAGS})
SET(BUILD_INFO_CXX_FLAGS ${CMAKE_CXX_FLAGS})
add_definitions(-D BUILD_INFO_C_FLAGS='"${BUILD_INFO_C_FLAGS}"')
add_definitions(-D BUILD_INFO_CXX_FLAGS='"${BUILD_INFO_CXX_FLAGS}"')


add_library(build_info STATIC
            build_info.cpp)

add_executable(test_build_info test_build_info.cpp)
target_link_libraries(test_build_info build_info)

include(CTest)
enable_testing()

add_custom_target(__alltests__)
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} -V)
add_dependencies(check __alltests__)

set(TEST_LIBRARIES pthread gmock_main gtest gmock)
set(TEST_WORKING_DIRECTORY ${IGREC_MAIN_SRC_DIR}/..)
set(TEST_COMMAND_ARGS --gtest_color=yes)

find_program(MEMORYCHECK_COMMAND valgrind)
if(MEMORYCHECK_COMMAND)
    set(MEMORYCHECK_COMMAND_OPTIONS "--trace-children=yes --leak-check=full")
    add_custom_target(memcheck COMMAND ${CMAKE_CTEST_COMMAND} -D NightlyMemCheck -V)
    add_dependencies(memcheck __alltests__)
endif()

function(make_test test)
  if (NOT TARGET ${test})
    add_executable(${test} ${ARGN})
    target_link_libraries(${test} ${TEST_LIBRARIES})
    get_property(BINDIR TARGET ${test} PROPERTY LOCATION)
    add_test(NAME ${test}
             WORKING_DIRECTORY ${TEST_WORKING_DIRECTORY}
             COMMAND ${test} ${TEST_COMMAND_ARGS})
  endif(NOT TARGET ${test})

  add_dependencies(__alltests__ ${test})
  set_target_properties(${test} PROPERTIES EXCLUDE_FROM_ALL 1)
endfunction(make_test)

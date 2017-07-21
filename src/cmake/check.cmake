include(CTest)

find_program(MEMORYCHECK_COMMAND valgrind)
if(NOT MEMORYCHECK_COMMAND)
  # add_dependencies(memcheck valgrind) It should be done further
  set(COMPILE_VALGRIND TRUE)
  get_target_property(MEMORYCHECK_COMMAND valgrind IMPORTED_LOCATION)
  message("valgrind not found, to be build at ${MEMORYCHECK_COMMAND}")
else()
  message("valgrind found at ${MEMORYCHECK_COMMAND}")
endif()

enable_testing()
add_custom_target(__alltests__)
add_custom_target(check COMMAND OMP_NUM_THREADS=1 ${CMAKE_CTEST_COMMAND} -V)
add_dependencies(check __alltests__)

set(TEST_LIBRARIES pthread gmock_main gtest gmock)
set(TEST_WORKING_DIRECTORY ${IGREC_MAIN_SRC_DIR}/..)
set(TEST_COMMAND_ARGS --gtest_color=yes)

set(MEMORYCHECK_COMMAND_OPTIONS "--trace-children=yes --leak-check=full")
add_custom_target(memcheck COMMAND OMP_NUM_THREADS=1 ${CMAKE_CTEST_COMMAND} -D NightlyMemCheck -V)
add_dependencies(memcheck __alltests__)
if(COMPILE_VALGRIND)
  add_dependencies(memcheck valgrind)
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

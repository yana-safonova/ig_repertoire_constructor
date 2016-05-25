include(CTest)
enable_testing()

add_custom_target(check_essential)
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} -V)

set(TEST_LIBRARIES pthread gtest_main gtest gmock)
set(TEST_WORKING_DIRECTORY ${SPADES_MAIN_SRC_DIR}/..)
set(TEST_COMMAND_ARGS --gtest_color=yes)

function(make_essential_test test)
  if (NOT TARGET ${test})
    add_executable(${test} ${ARGN})
    target_link_libraries(${test} ${TEST_LIBRARIES})
    get_property(BINDIR TARGET ${test} PROPERTY LOCATION)
    add_test(NAME ${test}
             WORKING_DIRECTORY ${TEST_WORKING_DIRECTORY}
             COMMAND ${test} ${TEST_COMMAND_ARGS})
  endif(NOT TARGET ${test})

  add_dependencies(check ${test})
  add_dependencies(check_essential ${test})
  add_custom_command(TARGET ${test} POST_BUILD COMMAND ${CMAKE_CTEST_COMMAND} -V -R ^${test}$$ || (${CMAKE_COMMAND} -E remove $<TARGET_FILE:${test}> && false))
  set_target_properties(${test} PROPERTIES EXCLUDE_FROM_ALL 1)
endfunction(make_essential_test)

function(make_test test)
  if (NOT TARGET ${test})
    add_executable(${test} ${ARGN})
    target_link_libraries(${test} ${TEST_LIBRARIES})
    get_property(BINDIR TARGET ${test} PROPERTY LOCATION)
    add_test(NAME ${test}
             WORKING_DIRECTORY ${TEST_WORKING_DIRECTORY}
             COMMAND ${test} ${TEST_COMMAND_ARGS})
  endif(NOT TARGET ${test})

  add_dependencies(check ${test})
  set_target_properties(${test} PROPERTIES EXCLUDE_FROM_ALL 1)
endfunction(make_test)

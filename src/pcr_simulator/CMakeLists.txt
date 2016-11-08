add_executable(simulate_barcoded
               main.cpp
               pcr_simulator.cpp
               pcr_simulator.hpp
               options.hpp)

target_link_libraries(simulate_barcoded input boost_program_options boost_filesystem boost_system ${COMMON_LIBRARIES})
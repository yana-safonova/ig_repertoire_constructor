function(enable_profiler)
    if (not ${GOOGLE_PERFTOOLS_FOUND})
        find_package(GooglePerfTools REQUIRED COMPONENTS tcmalloc profiler tcmalloc_and_profiler)
    endif()

    foreach(target ${ARGV})
        get_target_property(libs ${target} LINK_LIBRARIES)

        if ("${libs}" MATCHES "(profiler)|(tcmalloc_and_profiler)")
            # Do nothing
        else()
            if (";${libs};" MATCHES "tcmalloc")
                target_link_libraries(${target} "-Wl,--no-as-needed -ltcmalloc_and_profiler -Wl,--as-needed")
            else()
                target_link_libraries(${target} "-Wl,--no-as-needed -lprofiler -Wl,--as-needed")
            endif()
        endif()

        # get_target_property(libs ${target} LINK_LIBRARIES)
        # message(${libs})
    endforeach()
endfunction()

project(configs)
cmake_minimum_required(VERSION 2.8)

function(cpcfg overwrite)
    if (${overwrite})
        MESSAGE("Configs will be overwritten")
    endif()
    file(GLOB config_templates "*/*.template")
    foreach(config_template ${config_templates})
    # Remove extention
    STRING(REGEX REPLACE "\\.template$" "" config ${config_template})
    if(NOT EXISTS "${config}" OR ${overwrite})
        configure_file("${config_template}"
                        "${config}"
                        @ONLY)
    endif()
    endforeach()
endfunction(cpcfg)

if (NOT CMAKE_OVERWRITE_CONFIGS)
    SET(CMAKE_OVERWRITE_CONFIGS "FALSE")
endif()

cpcfg(${CMAKE_OVERWRITE_CONFIGS})

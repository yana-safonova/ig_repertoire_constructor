project(configs)
cmake_minimum_required(VERSION 2.8)

file(GLOB config_templates "*/*.template")
foreach(config_template ${config_templates})
  # Remove extention
  STRING(REGEX REPLACE "\\.template$" "" config ${config_template})
  if(NOT EXISTS "${config}")
    configure_file("${config_template}"
                    "${config}"
                    @ONLY)
  endif()
endforeach()

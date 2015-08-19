#pragma once


#include <string>


std::string join_cmd_line(size_t argc, char **argv) {
  std::string result = argv[0];
  for (size_t i = 1; i < argc; ++i) {
    result += " ";
    result += argv[i];
  }

  return result;
}

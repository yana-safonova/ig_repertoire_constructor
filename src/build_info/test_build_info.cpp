#include <build_info.hpp>

#include <iostream>

#define OUT_BUILD_INFO_FIELD(field) #field << ": " << build_info::field

int main() {
    std::cout << OUT_BUILD_INFO_FIELD(version) << std::endl;

    std::cout << OUT_BUILD_INFO_FIELD(git_hash) << std::endl;
    std::cout << OUT_BUILD_INFO_FIELD(git_hash7) << std::endl;
    std::cout << OUT_BUILD_INFO_FIELD(git_refspec) << std::endl;
    std::cout << OUT_BUILD_INFO_FIELD(git_describe) << std::endl;
    std::cout << OUT_BUILD_INFO_FIELD(git_tag) << std::endl;

    std::cout << OUT_BUILD_INFO_FIELD(cmake) << std::endl;
    std::cout << OUT_BUILD_INFO_FIELD(system) << std::endl;
    std::cout << OUT_BUILD_INFO_FIELD(cpu) << std::endl;

    std::cout << OUT_BUILD_INFO_FIELD(c_compiler) << std::endl;
    std::cout << OUT_BUILD_INFO_FIELD(c_compiler_id) << std::endl;
    std::cout << OUT_BUILD_INFO_FIELD(c_compiler_version) << std::endl;
    std::cout << OUT_BUILD_INFO_FIELD(c_flags) << std::endl;

    std::cout << OUT_BUILD_INFO_FIELD(cxx_compiler) << std::endl;
    std::cout << OUT_BUILD_INFO_FIELD(cxx_compiler_id) << std::endl;
    std::cout << OUT_BUILD_INFO_FIELD(cxx_compiler_version) << std::endl;
    std::cout << OUT_BUILD_INFO_FIELD(cxx_flags) << std::endl;

    std::cout << OUT_BUILD_INFO_FIELD(build_type) << std::endl;

    return 0;
}

// vim: ts=4:sw=4

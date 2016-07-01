#include <build_info.hpp>

#include <iostream>

int main() {
    std::cout << build_info::version << std::endl;

    std::cout << build_info::git_hash << std::endl;
    std::cout << build_info::git_hash7 << std::endl;
    std::cout << build_info::git_refspec << std::endl;
    std::cout << build_info::git_describe << std::endl;
    std::cout << build_info::git_tag << std::endl;

    std::cout << build_info::cmake << std::endl;
    std::cout << build_info::system << std::endl;
    std::cout << build_info::cpu << std::endl;

    std::cout << build_info::c_compiler_id << std::endl;
    std::cout << build_info::c_compiler << std::endl;
    std::cout << build_info::cxx_compiler_id << std::endl;
    std::cout << build_info::cxx_compiler << std::endl;

    std::cout << build_info::build_type << std::endl;
    std::cout << build_info::c_flags << std::endl;
    std::cout << build_info::cxx_flags << std::endl;

    return 0;
}

// vim: ts=4:sw=4

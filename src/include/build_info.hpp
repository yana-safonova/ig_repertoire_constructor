#pragma once

namespace build_info {
    extern const char *version;

    extern const char *git_hash;
    extern const char *git_hash7;
    extern const char *git_refspec;
    extern const char *git_describe;
    extern const char *git_tag;

    extern const char *cmake;
    extern const char *system;
    extern const char *cpu;

    extern const char *c_compiler_id;
    extern const char *cxx_compiler_id;
    extern const char *c_compiler;
    extern const char *cxx_compiler;

    extern const char *build_type;
    extern const char *c_flags;
    extern const char *cxx_flags;
} // namespace build_info

// vim: ts=4:sw=4

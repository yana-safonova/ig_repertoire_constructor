//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once
#include "stacktrace.hpp"

#include "boost/current_function.hpp"
#include <sstream>
#include <iostream>
#include <cstdio>

#ifndef NVERIFY
#define VERIFY(expr)                                                        \
    do {                                                                    \
        if (!(expr)) {                                                      \
            std::stringstream ss;                                           \
            print_stacktrace();                                             \
            ss << "Verification of expression '" << #expr <<                \
                    "' failed in function '" <<  BOOST_CURRENT_FUNCTION <<  \
                    "'. In file '" << __FILE__ <<                           \
                    "' on line " << __LINE__ << ".";                        \
            std::cout << ss.str() << std::endl;                             \
            std::cerr << ss.str() << std::endl;                             \
            fflush(stdout);                                                 \
            fflush(stderr);                                                 \
            abort();                                                        \
        }                                                                   \
    } while(0);

#define VERIFY_MSG(expr, msg)                                               \
    do {                                                                    \
        if (!(expr)) {                                                      \
            std::stringstream ss;                                           \
            print_stacktrace();                                             \
            ss << "Verification of expression '" << #expr <<                \
                    "' failed in function '" <<  BOOST_CURRENT_FUNCTION <<  \
                    "'. In file '" << __FILE__ <<                           \
                    "' on line " << __LINE__ <<                             \
                    ". Message '" << msg << "'.";                           \
            std::cout << ss.str() << std::endl;                             \
            std::cerr << ss.str() << std::endl;                             \
            fflush(stdout);                                                 \
            fflush(stderr);                                                 \
            abort();                                                        \
        }                                                                   \
    } while(0);
#else
#define VERIFY(expr) ((void) 0);
#define VERIFY_MSG(expr, msg) ((void) 0);
#endif

namespace verify_details {

template<typename T>
void __FormVerifyMessage(std::ostringstream &ss, T &&t) {
    ss << t;
}

template<typename T, typename... Args>
void __FormVerifyMessage(std::ostringstream &ss, T &&t, Args &&... args) {
    ss << t << " ";
    __FormVerifyMessage(ss, std::forward<Args>(args)...);
}

} // End namespace verify_details

template<typename... Args>
std::string FormVerifyMessage(Args &&... args) {
    std::ostringstream ss;
    verify_details::__FormVerifyMessage(ss, std::forward<Args>(args)...);
    return ss.str();
}

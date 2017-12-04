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
            std::stringstream ss_verify;                                    \
            print_stacktrace();                                             \
            ss_verify << "Verification of expression '" << #expr <<         \
                    "' failed in function '" <<  BOOST_CURRENT_FUNCTION <<  \
                    "'. In file '" << __FILE__ <<                           \
                    "' on line " << __LINE__ << ".";                        \
            std::cout << ss_verify.str() << std::endl;                      \
            std::cerr << ss_verify.str() << std::endl;                      \
            fflush(stdout);                                                 \
            fflush(stderr);                                                 \
            abort();                                                        \
        }                                                                   \
    } while(0);

#define VERIFY_MSG(expr, msg)                                               \
    do {                                                                    \
        if (!(expr)) {                                                      \
            std::stringstream ss_verify;                                    \
            print_stacktrace();                                             \
            ss_verify << "Verification of expression '" << #expr <<         \
                    "' failed in function '" <<  BOOST_CURRENT_FUNCTION <<  \
                    "'. In file '" << __FILE__ <<                           \
                    "' on line " << __LINE__ <<                             \
                    ". Message '" << msg << "'.";                           \
            std::cout << ss_verify.str() << std::endl;                      \
            std::cerr << ss_verify.str() << std::endl;                      \
            fflush(stdout);                                                 \
            fflush(stderr);                                                 \
            abort();                                                        \
        }                                                                   \
    } while(0);
#else
#define VERIFY(expr) ((void) 0);
#define VERIFY_MSG(expr, msg) ((void) 0);
#endif

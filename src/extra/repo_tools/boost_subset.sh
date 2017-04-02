#!/bin/bash

REQUIRED_BOOST_LIBRARIES=" regex filesystem random system iostreams program_options "
REQUIRED_BOOST_LIBRARIES+=" format algorithm  property_tree tokenizer bind function optional lexical_cast noncopyable foreach iterator "
REQUIRED_BOOST_LIBRARIES+=" utility type_traits current_function pending "
REQUIRED_BOOST_LIBRARIES+=" pool "
REQUIRED_BOOST_LIBRARIES+=" multi_array "
# REQUIRED_BOOST_LIBRARIES+="atomic "  # needed for thread

EXCLUDED_BINARIES=" config serialization thread date_time exception smart_ptr test timer "

BOOST_VERSION="1_63"
# BOOST_VERSION="1_59"

OUTPUT_DIR="out_${BOOST_VERSION}_0"
BOOST_DIR="boost_${BOOST_VERSION}_0"

# INCLUDE_DIR=$(realpath ${OUTPUT_DIR})
INCLUDE_DIR="\${EXT_DIR}/include"

mkdir -p "${OUTPUT_DIR}"

bcp --boost="${BOOST_DIR}" --unix-lines ${REQUIRED_BOOST_LIBRARIES} "${OUTPUT_DIR}"
bcp --boost="${BOOST_DIR}" --list ${REQUIRED_BOOST_LIBRARIES} > "${OUTPUT_DIR}/filelist.txt"
cd ${OUTPUT_DIR}
rm -fr doc
rm Jamroot boost.png rst.css boost.css

SRC_DIR="boost_src"
mkdir -p "${SRC_DIR}"
echo "# -*- cmake -*-

cmake_minimum_required(VERSION 2.8)

add_definitions(-w)

include_directories(${INCLUDE_DIR})

" > "${SRC_DIR}/CMakeLists.txt"

for library_source in libs/*
do
    library=${library_source##*/}
    echo "Processing ${library} in ${library_source}..."
    if [ ! -d "libs/${library}/src" ]; then
        echo "Not a compileable lib!"
        continue
    fi
    if (echo "${EXCLUDED_BINARIES}"  | fgrep -q -w "${library}") then
        echo "${library} excluded!"
        continue
    fi
    rm -fr "${SRC_DIR}/${library}"
    cp -R "libs/$library/src/" "${SRC_DIR}/${library}"
    echo "copied libs/$library/src/ -> ${SRC_DIR}/${library}"

    echo "add_subdirectory($library)" >> "${SRC_DIR}/CMakeLists.txt"
    echo "cmake_minimum_required(VERSION 2.8)

project(boost_$library CXX)

file(GLOB_RECURSE boost_${library}_source_files \"*.cpp\")
add_library(boost_${library} STATIC
            \${boost_${library}_source_files})" > "${SRC_DIR}/${library}/CMakeLists.txt"
done

rm -fr "${SRC_DIR}/thread/win32"
rm "${SRC_DIR}/test/cpp_main.cpp"
mv "${SRC_DIR}/thread/pthread/once_atomic.cpp" "${SRC_DIR}/thread/pthread/once_atomic.inc"
sed -i s/once_atomic.cpp/once_atomic.inc/ ${SRC_DIR}/thread/pthread/once.cpp

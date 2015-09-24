project(bzip2 C)

include_directories(${EXT_DIR}/include)

file(GLOB bzip2_SRC
     "*.c"
)

add_library(bz2 STATIC
            ${bzip2_SRC})


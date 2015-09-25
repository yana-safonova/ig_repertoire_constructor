project(zlib C)

include_directories(${EXT_DIR}/include)

file(GLOB zlib_SRC
     "*.c"
)

add_library(z STATIC
            ${zlib_SRC})

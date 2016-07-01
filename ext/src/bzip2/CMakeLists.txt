project(bzip2 C)

include_directories(${EXT_DIR}/include)

add_library(bz2 STATIC
            blocksort.c
            huffman.c
            crctable.c
            randtable.c
            compress.c
            decompress.c
            bzlib.c)

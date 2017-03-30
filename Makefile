# Default build type
build_type?="RelWithAsserts"

# Default install prefix
prefix?="/usr/local"

.PHONY: clean cleanup cmake

all: igrec

cmake:
	mkdir -p build/release
	cd build/release && cmake ../.. -DCMAKE_BUILD_TYPE="${build_type}" -DCMAKE_INSTALL_PREFIX=${prefix} -Wno-dev

igrec: cmake
	$(MAKE) -C build/release all

pack: igrec
	cd build/release && cpack

deb: igrec
	cd build/release && cpack -G DEB

install: igrec
	cd build/release && cmake -P cmake_install.cmake

rig: cmake
	$(MAKE) -C build/release/ig_repertoire_constructor ig_repertoire_constructor

dsf: cmake
	$(MAKE) -C build/release/dense_sgraph_finder dense_sgraph_finder

check: cmake
	$(MAKE) -C build/release check

memcheck: cmake
	$(MAKE) -C build/release memcheck

rnd: cmake
	$(MAKE) -C build/release memcheck

vjf: cmake
	$(MAKE) -C build/release/vj_finder

cdr: cmake
	$(MAKE) -C build/release/cdr_labeler

umi: cmake
	$(MAKE) -C build/release/umi_experiments

clean:
	-rm -fr build/release

cleanup:
	rm *.pyc
	rm -rf igrec_test
	rm -rf ms_analyzer_test
	rm *~

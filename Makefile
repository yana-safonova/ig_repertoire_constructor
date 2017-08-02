# Default build type
build_type?="RelWithAsserts"

# Default install prefix
prefix?="/usr/local"

# Not-verify (disabled)
nverify?=""

# Build with gpertools profiler (disabled)
gperf?=""

.PHONY: clean clean_tests cmake all pack

all: igrec

cpcfg:
	mkdir -p build/tmp
	cd build/tmp && cmake ../../configs -DCMAKE_OVERWRITE_CONFIGS=true
	rm -r build/tmp

cmake:
	mkdir -p build/release
	cd build/release && cmake ../.. -DCMAKE_BUILD_TYPE="${build_type}" -DCMAKE_INSTALL_PREFIX=${prefix} \
		-Wno-dev \
		-DCMAKE_NVERIFY=${nverify} \
		-DCMAKE_GOOGLE_PROFILER=${gperf}

igrec: cmake
	$(MAKE) -C build/release all

rpm: igrec
	cd build/release && cpack -G RPM

deb: igrec
	cd build/release && cpack -G DEB

tgz: igrec
	cd build/release && cpack -G TGZ

install: igrec
	cd build/release && cmake -P cmake_install.cmake

dsf: cmake
	$(MAKE) -C build/release dense_sgraph_finder

shm_kmer_matrix: cmake
	$(MAKE) -C build/release/shm_kmer_matrix_estimator

check: cmake
	$(MAKE) -C build/release check

memcheck: cmake
	$(MAKE) -C build/release memcheck

rnd: cmake
	$(MAKE) -C build/release rnd

vjf: cmake
	$(MAKE) -C build/release vj_finder

cdr: cmake
	$(MAKE) -C build/release cdr_labeler

umi: cmake
	$(MAKE) -C build/release umi_correction_stats umi_graph umi_naive umi_to_fastq

clean:
	-rm *.pyc
	-rm py/*.pyc
	-rm -r build

ant: cmake
	$(MAKE) -C build/release antevolo

igs: cmake
	$(MAKE) -C build/release ig_simulator

clean_tests:
	-rm *.pyc
	-rm -r igrec_test
	-rm -r ms_analyzer_test
	-rm -r ig_simulator_test
	-rm *~

all: igrec

igrec: build/release/Makefile
	$(MAKE) -C build/release all

build/release/Makefile: prepare_cfg
	./prepare_cfg

# Default install prefix
prefix?="/opt/"

install: igrec
	cd build/release && cmake -DCMAKE_INSTALL_PREFIX=${prefix} -P cmake_install.cmake

rig:
	$(MAKE) -C build/release/ig_repertoire_constructor ig_repertoire_constructor

dig:
	$(MAKE) -C build/debug/ig_repertoire_constructor ig_repertoire_constructor

dsf:
	$(MAKE) -C build/release/dense_sgraph_finder dense_sgraph_finder

metis:
	$(MAKE) -C build/release/ext_tools/metis-5.1.0/ metis

check: all
	$(MAKE) -C build/release check

check_essential: all
	$(MAKE) -C build/release check_essential

all:
	$(MAKE) -C build/release all

rig:
	$(MAKE) -C build/release/ig_repertoire_constructor ig_repertoire_constructor

dig:
	$(MAKE) -C build/debug/ig_repertoire_constructor ig_repertoire_constructor

dsf:
	$(MAKE) -C build/release/dense_sgraph_finder dense_sgraph_finder

metis:
	$(MAKE) -C build/release/ext_tools/metis-5.1.0/ metis

shm:
	$(MAKE) -C build/release/shm_kmer_model/ shm_kmer_model

check: all
	$(MAKE) -C build/release check

check_essential: all
	$(MAKE) -C build/release check_essential

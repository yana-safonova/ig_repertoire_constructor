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

shm_kmer_matrix:
	$(MAKE) -C build/release/shm_kmer_matrix_estimator

check:
	$(MAKE) -C build/release check

memcheck:
	$(MAKE) -C build/release memcheck

vjf:
	$(MAKE) -C build/release/vj_finder

cdr:
	$(MAKE) -C build/release/cdr_labeler

ant:
	$(MAKE) -C build/release/antevolo

umi:
	$(MAKE) -C build/release/umi_experiments

all:
	$(MAKE) -C build/release all
	$(MAKE) -C build/debug all

rig:
	$(MAKE) -C build/release/ig_repertoire_constructor ig_repertoire_constructor

dig:
	$(MAKE) -C build/debug/ig_repertoire_constructor ig_repertoire_constructor

igtools:
	g++ src/ig_tools/paired_read_merger/main.cpp -o build/release/bin/paired_read_merger
	g++ src/ig_tools/fastq_to_fasta/fastq_to_fasta.cpp -o build/release/bin/fastq_to_fasta
	g++ -std=c++11 src/ig_tools/merged_reads_stats_calculator/main.cpp -o build/release/bin/compute_merged_reads_stats

dsf:
	$(MAKE) -C build/release/dense_sgraph_finder dense_sgraph_finder

metis:
	$(MAKE) -C build/release/ext_tools/metis-5.1.0/ metis

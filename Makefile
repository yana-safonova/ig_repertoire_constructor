all:
	$(MAKE) -C build/release all
	$(MAKE) -C build/debug all
	g++ src/ig_tools/paired_read_merger/main.cpp -o bin/ig_tools/paired_read_merger
	g++ src/ig_tools/fastq_to_fasta/fastq_to_fasta.cpp -o bin/ig_tools/fastq_to_fasta
	g++ -std=c++11 src/ig_tools/merged_reads_stats_calculator/main.cpp -o bin/ig_tools/compute_merged_reads_stats

rig:
	$(MAKE) -C build/release/ig_repertoire_constructor ig_repertoire_constructor

dig:
	$(MAKE) -C build/debug/ig_repertoire_constructor ig_repertoire_constructor

igtools:
	g++ src/ig_tools/paired_read_merger/main.cpp -o bin/ig_tools/paired_read_merger
	g++ src/ig_tools/fastq_to_fasta/fastq_to_fasta.cpp -o bin/ig_tools/fastq_to_fasta
	g++ -std=c++11 src/ig_tools/merged_reads_stats_calculator/main.cpp -o bin/ig_tools/compute_merged_reads_stats

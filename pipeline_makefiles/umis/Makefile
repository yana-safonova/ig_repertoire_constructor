include ../Makefile_vars

VJF = ../vj_finder

all:
	$(MAKE) -C $(VJF)
	$(MAKE) cleaned_reads_umi.fastq cleaned_reads_umi_compressed.fastq cleaned_reads_umi_compressed.graph

cleaned_reads_umi.fastq cleaned_reads_umi_compressed.fastq cleaned_reads_umi_compressed.graph: $(VJF)/cleaned_reads.fa !! $(UMI_SRC)/stats/print_umi_graph_stats.py
	$(UMI_SRC)/stats/print_umi_graph_stats.py \
	--input $(VJF)/cleaned_reads.fa \
	--output stats \
	--tmp . \
	--tau %UMI_GRAPH_TAU \
	--clean

clean:
	$(MAKE) -C $(VJF) clean
	$(RM) cleaned_reads_umi.fastq cleaned_reads_umi_compressed.fastq cleaned_reads_umi_compressed.graph

include ../Makefile_vars

VJF = ../vj_finder
UMIS = ../umis

all:
	$(MAKE) -C $(VJF)
	$(MAKE) -C $(UMIS)
	$(MAKE) $(INTER_REP)

$(INTER_REP): $(BIN)/cluster_reads $(VJF)/cleaned_reads.fa $(UMIS)/cleaned_reads_umi.fastq $(UMIS)/cleaned_reads_umi_compressed.fastq $(UMIS)/cleaned_reads_umi_compressed.graph
	$(BIN)/cluster_reads \
	--reads $(VJF)/cleaned_reads.fa \
	--umi_uncompressed $(UMIS)/cleaned_reads_umi.fastq \
	--umi_compressed $(UMIS)/cleaned_reads_umi_compressed.fastq \
	--graph $(UMIS)/cleaned_reads_umi_compressed.graph \
	--output . \
	--save-clusters false

clean:
	$(MAKE) -C $(VJF) clean
	$(MAKE) -C $(UMIS) clean
	$(RM) intermediate_repertoire.*
	$(RM) -r clusters_by_umis

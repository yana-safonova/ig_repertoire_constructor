include ../Makefile_vars

VJF = ../vj_finder
UMIS = ../umis

all:
	$(MAKE) -C $(VJF)
	$(MAKE) -C $(UMIS)
	$(MAKE) $(INTER_REP)

$(INTER_REP): $(BIN)/cluster_reads $(VJF)/cleaned_reads.fa $(UMIS)/cleaned_reads_umi.fastq $(UMIS)/cleaned_reads_umi_compressed.fastq $(UMIS)/cleaned_reads_umi_compressed.graph
	$(BIN)/cluster_reads \
	-r $(VJF)/cleaned_reads.fa \
	-u $(UMIS)/cleaned_reads_umi.fastq \
	-c $(UMIS)/cleaned_reads_umi_compressed.fastq \
	-g $(UMIS)/cleaned_reads_umi_compressed.graph \
	-o . \
	-s true

clean:
	$(MAKE) -C $(VJF) clean
	$(MAKE) -C $(UMIS) clean
	$(RM) intermediate_repertoire.*
	$(RM) -r clusters_by_umis

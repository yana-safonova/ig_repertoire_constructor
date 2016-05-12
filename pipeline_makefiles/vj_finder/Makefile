include ../Makefile_vars

MERGED = ../merged_reads

all:
	$(MAKE) -C $(MERGED)
	$(MAKE) cleaned_reads.fa

cleaned_reads.fa: $(BIN)/ig_kplus_vj_finder $(MERGED)/merged_reads.fastq
	$(BIN)/ig_kplus_vj_finder \
	--input $(MERGED)/merged_reads.fastq \
	-o . \
	--db-directory $(IGREC)/data/germline \
	-t %THREADS \
	--loci %LOCI \
	--organism %ORGANISM

clean:
	$(MAKE) -C $(MERGED) clean
	$(RM) alignment_info.csv cleaned_reads.fa filtered_reads.fa

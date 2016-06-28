include ../Makefile_vars

MERGED = ../merged_reads

all:
	$(MAKE) -C $(MERGED)
	$(MAKE) cleaned_reads.fa

cleaned_reads.fa: $(MERGED)/merged_reads.fastq $(BIN)/ig_kplus_vj_finder
	$(BIN)/ig_kplus_vj_finder \
	--input-file $(MERGED)/merged_reads.fastq \
	--output . \
	--db-directory $(IGREC)/data/germline \
	--threads %THREADS \
	--loci %LOCI \
	--organism %ORGANISM

clean:
	$(MAKE) -C $(MERGED) clean
	$(RM) alignment_info.csv cleaned_reads.fa filtered_reads.fa

include ../Makefile_vars

MERGED = ../merged_reads

all:
	$(MAKE) -C $(MERGED)
	$(MAKE) cleaned_reads.fa

cleaned_reads.fa: $(MERGED)/merged_reads.fastq !! $(BIN)/$(VJ_FINDER)
	cd $(IGREC); \
	$(BIN)/$(VJ_FINDER) \
	--input-file $(CURDIR)/$(MERGED)/merged_reads.fastq \
	--output-dir $(CURDIR)/ \
	--db-directory $(IGREC)/data/germline \
	--threads %THREADS \
	--loci %LOCI \
	--organism %ORGANISM

clean:
	$(MAKE) -C $(MERGED) clean
	$(RM) alignment_info.csv cleaned_reads.fa filtered_reads.fa

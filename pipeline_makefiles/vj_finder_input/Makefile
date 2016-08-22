include ../Makefile_vars

all:
	$(MAKE) cleaned_reads.fa

cleaned_reads.fa: %INPUT !! $(BIN)/$(VJ_FINDER)
	cd $(IGREC); \
	$(BIN)/$(VJ_FINDER) \
	--input-file %INPUT \
	--output-dir $(CURDIR)/ \
	--db-directory $(IGREC)/data/germline \
	--threads %THREADS \
	--loci %LOCI \
	--organism %ORGANISM

clean:
	$(RM) alignment_info.csv cleaned_reads.fa filtered_reads.fa

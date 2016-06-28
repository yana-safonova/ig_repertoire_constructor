include ../Makefile_vars

all:
	$(MAKE) cleaned_reads.fa

cleaned_reads.fa: %INPUT $(BIN)/ig_kplus_vj_finder
	$(BIN)/ig_kplus_vj_finder \
	--input-file %INPUT \
	--output . \
	--db-directory $(IGREC)/data/germline \
	--threads %THREADS \
	--loci %LOCI \
	--organism %ORGANISM

clean:
	$(RM) alignment_info.csv cleaned_reads.fa filtered_reads.fa

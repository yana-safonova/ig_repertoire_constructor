include ../Makefile_vars

all:
	$(MAKE) cleaned_reads.fa

cleaned_reads.fa: $(BIN)/ig_kplus_vj_finder %INPUT
	$(BIN)/ig_kplus_vj_finder \
	--input %INPUT \
	-o . \
	--db-directory $(IGREC)/data/germline \
	-t %THREADS \
	--loci %LOCI \
	--organism %ORGANISM

clean:
	$(RM) alignment_info.csv cleaned_reads.fa filtered_reads.fa

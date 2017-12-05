include ../Makefile_vars

VJ_FINDER_COMP = $(CURDIR)/../vj_finder_compressed
COMPRESSED_INPUT := $(VJ_FINDER_COMP)/compressed%EXT_INPUT

all:
	$(MAKE) cleaned_reads.fa

cleaned_reads.fa: %INPUT !! $(BIN)/$(VJ_FINDER)
	mkdir -p $(VJ_FINDER_COMP)

	$(BIN)/compress_ion \
	--input %INPUT \
	--compressed $(COMPRESSED_INPUT)

	cd $(IGREC); \
	$(BIN)/$(VJ_FINDER) \
	--input-file $(COMPRESSED_INPUT) \
	--output-dir $(VJ_FINDER_COMP) \
	--db-directory $(IGREC)/data/germline \
	--threads %THREADS \
	--loci %LOCI \
	--organism human_ion \
	--left-uncoverage-limit 30 \
	--"min-vsegment-length" 150 \
	--min-jsegment-length 20 \
	--min-len 180 \
	--fix-left 0 \
	--fix-right 0

	$(BIN)/uncompress_vjf_results \
	--cleaned $(VJ_FINDER_COMP)/cleaned_reads.fa \
	--align $(VJ_FINDER_COMP)/alignment_info.csv \
	--reads %INPUT \
	--v $(GERMLINE_ION)/IGHV.fa \
	--j $(GERMLINE_ION)/IGHJ.fa \
	--output cleaned_reads.fa

clean:
	$(RM) alignment_info.csv cleaned_reads.fa filtered_reads.fa

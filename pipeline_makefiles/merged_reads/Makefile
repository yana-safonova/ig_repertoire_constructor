include ../Makefile_vars

all:
	$(MAKE) merged_reads.fastq

merged_reads.fastq: %LEFT_INPUT %RIGHT_INPUT !! $(BIN)/paired_read_merger
	$(BIN)/paired_read_merger \
	%LEFT_INPUT \
	%RIGHT_INPUT \
	merged_reads.fastq

clean:
	$(RM) merged_reads.fastq

include ../Makefile_vars

all:
	$(MAKE) merged_reads.fastq

merged_reads.fastq: $(BIN)/paired_read_merger
	realpath ../age_ig_s3_R1.fastq
	$(BIN)/paired_read_merger \
	../age_ig_s3_R1.fastq \
	../age_ig_s3_R2.fastq \
	merged_reads.fastq

clean:
	$(RM) merged_reads.fastq

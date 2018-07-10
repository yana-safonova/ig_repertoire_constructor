include ../Makefile_vars

all:
	$(MAKE) cleaned_reads.fa

cleaned_reads.fa: %INPUT
	cp %INPUT $(CURDIR)/cleaned_reads.fa

clean:
	$(RM) alignment_info.csv cleaned_reads.fa filtered_reads.fa

include ../Makefile_vars

all:
	$(MAKE) -C ../compilation
	$(MAKE) -C ../umi_clustering
	$(MAKE) final_repertoire.fa

final_repertoire.fa: $(IGREC)/igrec.py ../umi_clustering/$(INTER_REP)
	$(IGREC)/igrec.py \
	-s ../umi_clustering/$(INTER_REP) \
	--loci %LOCI \
	--organism %ORGANISM \
	--tau %IGREC_TAU \
	--output . \
	--threads %THREADS

clean:
	$(MAKE) -C ../umi_clustering clean

clean_all:
	$(MAKE) clean
	$(MAKE) -C ../compilation clean

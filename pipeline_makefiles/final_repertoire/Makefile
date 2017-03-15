include ../Makefile_vars

CONSENSUS = ../ig_consensus_finder
UMI_CLUSTERING = ../umi_clustering

all:
	$(MAKE) -C ../compilation
	$(MAKE) -C ../ig_consensus_finder
	$(MAKE) final_repertoire.fa final_repertoire.rcm

final_repertoire.fa final_repertoire.rcm: $(CONSENSUS)/final_repertoire_uncompressed.fa $(CONSENSUS)/final_repertoire_uncompressed.rcm !! $(BIN)/ig_trie_compressor $(PY_SRC)/adhoc_merge_rcm_with_id_list.py $(PY_SRC)/fix_cluster_ids.py
	$(BIN)/ig_trie_compressor \
	--input-file $(CONSENSUS)/final_repertoire_uncompressed.fa \
	--idmap id_map.tmp \
	--output-file tmp.fa

	python $(PY_SRC)/adhoc_merge_rcm_with_id_list.py \
	$(CONSENSUS)/final_repertoire_uncompressed.rcm \
	id_map.tmp \
	$(CONSENSUS)/final_repertoire_uncompressed.fa \
	tmp.fa \
	final_repertoire.rcm

	python $(PY_SRC)/fix_cluster_ids.py \
	tmp.fa \
	final_repertoire.rcm \
	final_repertoire.fa

	$(BIN)/report_umi_abundance \
	--final-rcm final_repertoire.rcm \
	--inter-rcm $(UMI_CLUSTERING)/intermediate_repertoire_close_umi.rcm \
	--repertoire final_repertoire.fa \
	--output final_repertoire_umi.fa.gz

clean:
	$(MAKE) -C $(CONSENSUS) clean

clean_all:
	$(MAKE) clean
	$(MAKE) -C ../compilation clean

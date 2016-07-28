include ../Makefile_vars

UMI_CLUSTERING = ../umi_clustering

all:
	$(MAKE) -C $(UMI_CLUSTERING)
	$(MAKE) $(COMPRESSED_INTER_REP)

$(COMPRESSED_INTER_REP): $(UMI_CLUSTERING)/$(INTER_REP) $(UMI_CLUSTERING)/$(INTER_REP_RCM) !! $(BIN)/ig_trie_compressor $(PY_SRC)/adhoc_merge_rcm_with_id_list.py
	$(BIN)/ig_trie_compressor \
	--input-file $(UMI_CLUSTERING)/$(INTER_REP) \
	--idmap id_map.tmp \
	--output-file $(COMPRESSED_INTER_REP)

	python $(PY_SRC)/adhoc_merge_rcm_with_id_list.py \
	$(UMI_CLUSTERING)/$(INTER_REP_RCM) \
	id_map.tmp \
	$(UMI_CLUSTERING)/$(INTER_REP) \
	$(COMPRESSED_INTER_REP) \
	$(COMPRESSED_INTER_REP_RCM)

clean:
	$(RM) $(COMPRESSED_INTER_REP)
	$(RM) $(COMPRESSED_INTER_REP_RCM)
	$(MAKE) -C $(UMI_CLUSTERING) clean

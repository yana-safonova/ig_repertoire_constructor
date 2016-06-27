include ../Makefile_vars

UMI_CLUSTERING = ../umi_clustering

all:
	$(MAKE) -C $(UMI_CLUSTERING)
	$(MAKE) $(COMPRESSED_INTER_REP)

$(COMPRESSED_INTER_REP): $(PY_SRC)/ig_compress_equal_clusters.py $(BIN)/ig_trie_compressor $(UMI_CLUSTERING)/$(INTER_REP) $(UMI_CLUSTERING)/$(INTER_REP_RCM)
	$(BIN)/ig_trie_compressor \
	--input-file $(UMI_CLUSTERING)/$(INTER_REP) \
	--idmap $(COMPRESSED_INTER_REP_RCM) \
	--output-file $(COMPRESSED_INTER_REP)
#	python $(PY_SRC)/ig_compress_equal_clusters.py \
#	$(UMI_CLUSTERING)/$(INTER_REP) \
#	$(COMPRESSED_INTER_REP) \
#	--tmp-fa-file tmp.fa \
#	--tmp-map-file tmp.rcm \
#	--rcm $(UMI_CLUSTERING)/$(INTER_REP_RCM) \
#	--output-rcm $(COMPRESSED_INTER_REP_RCM)

clean:
	$(RM) $(COMPRESSED_INTER_REP)
	$(RM) $(COMPRESSED_INTER_REP_RCM)
	$(MAKE) -C $(UMI_CLUSTERING) clean

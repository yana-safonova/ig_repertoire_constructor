include ../Makefile_vars

COMPRESSED_INTER_PATH = ../intermediate_ig_trie_compressor
VJF = ../vj_finder
DSF_PATH = ../dense_subgraph_finder

all:
	$(MAKE) -C $(DSF_PATH)
	$(MAKE) final_repertoire_uncompressed.fa final_repertoire_uncompressed.rcm

final_repertoire_uncompressed.fa final_repertoire_uncompressed.rcm: $(VJF)/cleaned_reads.fa $(COMPRESSED_INTER_PATH)/$(COMPRESSED_INTER_REP_RCM) $(DSF_PATH)/$(DENSE_SUBGRAPHS) !! $(BIN)/ig_component_splitter $(PY_SRC)/adhoc_merge_rcm_with_id_list.py
	python $(PY_SRC)/adhoc_merge_rcm_with_id_list.py \
	$(COMPRESSED_INTER_PATH)/$(COMPRESSED_INTER_REP_RCM) \
	$(DSF_PATH)/$(DENSE_SUBGRAPHS) \
	$(COMPRESSED_INTER_PATH)/$(COMPRESSED_INTER_REP) \
	- \
	final_repertoire_uncompressed.rcm

	$(BIN)/ig_component_splitter \
	--input-file $(VJF)/cleaned_reads.fa \
	--rcm-file final_repertoire_uncompressed.rcm \
	--output-file final_repertoire_uncompressed.fa \
	--output-rcm-file final_repertoire_uncompressed_after_split.rcm \
	--allow-unassigned true

clean:
	$(RM) final_repertoire_uncompressed.fa final_repertoire_uncompressed.rcm
	$(MAKE) -C $(IG_GRAPH_CONSTRUCTOR) clean

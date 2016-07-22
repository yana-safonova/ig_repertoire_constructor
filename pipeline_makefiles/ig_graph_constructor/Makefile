include ../Makefile_vars

COMPRESSED_INTER_PATH = ../intermediate_ig_trie_compressor

all:
	$(MAKE) -C $(COMPRESSED_INTER_PATH)
	$(MAKE) $(SW_GRAPH)

$(SW_GRAPH): $(COMPRESSED_INTER_PATH)/$(COMPRESSED_INTER_REP) !! $(BIN)/ig_swgraph_construct
	$(BIN)/ig_swgraph_construct \
	--input-file $(COMPRESSED_INTER_PATH)/$(COMPRESSED_INTER_REP) \
	--output-file $(SW_GRAPH) \
	--threads %THREADS \
	--tau %IGREC_TAU \
	--export-abundances

clean:
	$(RM) $(SW_GRAPH)
	$(MAKE) -C $(COMPRESSED_INTER_PATH) clean

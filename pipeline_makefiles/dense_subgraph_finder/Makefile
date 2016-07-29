include ../Makefile_vars

IG_GRAPH_CONSTRUCTOR = ../ig_graph_constructor

all:
	$(MAKE) -C $(IG_GRAPH_CONSTRUCTOR)
	$(MAKE) $(DENSE_SUBGRAPHS)

$(DENSE_SUBGRAPHS): $(IG_GRAPH_CONSTRUCTOR)/$(SW_GRAPH) !! $(IGREC)/dense_subgraph_finder.py
	$(IGREC)/dense_subgraph_finder.py \
	--graph $(IG_GRAPH_CONSTRUCTOR)/$(SW_GRAPH) \
	--output . \
	--threads %THREADS \
	--min-snode-size %MIN_SUPER_NODE_SIZE \
	--min-fillin %MIN_FILLIN \
	--clean-output-dir False

clean:
	$(RM) $(DENSE_SUBGRAPHS)
	$(MAKE) -C $(IG_GRAPH_CONSTRUCTOR) clean

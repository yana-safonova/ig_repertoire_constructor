#!/usr/bin/env python2

import os
import sys

home_directory = os.path.abspath(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))  # + '/'
py_src = os.path.join(home_directory, "py/pipeline/")

sys.path.append(py_src)
sys.path.append(home_directory)
# import support

pprof_bin = home_directory + "/build/release/thirdparty/bin/pprof"

if __name__ == "__main__":
    dir = sys.argv[1]
    binaries = {"consensus_constructor": "ig_component_splitter",
                "vjf": "vj_finder",
                "trie_compressor": "ig_trie_compressor",
                "graph_constructor": "ig_swgraph_construct"}
    for tool in ["consensus_constructor", "vjf", "trie_compressor", "graph_constructor"]:
        ret = os.system("%(pprof_bin)s --pdf %(home_dir)s/build/release/bin/%(bin)s %(dir)s/%(tool)s_prof.out > %(dir)s/%(tool)s.pdf" %
                        {"pprof_bin": pprof_bin,
                         "home_dir": home_directory,
                         "dir": dir,
                         "tool": tool,
                         "bin": binaries[tool]})

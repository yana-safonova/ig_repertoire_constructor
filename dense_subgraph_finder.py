#!/usr/bin/env python

import os
import sys
import init

def main():
    # just draft of main fuction
    command_line = init.PathToBins.run_dense_sgraph_finder + " configs/dense_subgraph_finder/config.info"
    error_code = os.system(command_line)
    if error_code != 0:
        print "ERROR: Dense sgraph finder finished abnormally"
        sys.exit(1)

if __name__ == '__main__':
    main() 

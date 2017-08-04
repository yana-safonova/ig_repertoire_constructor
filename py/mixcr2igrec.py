#!/usr/bin/env python2

from simulate import convert_mixcr_output_to_igrec

import sys

if __name__ == "__main__":
    convert_mixcr_output_to_igrec(sys.argv[1], sys.argv[2])

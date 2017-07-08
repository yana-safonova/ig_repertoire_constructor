#!/usr/bin/env python2

import os
import sys

input_dir = sys.argv[1]
files = os.listdir(input_dir)
for f in files:
    os.system('dot -Tpdf ' + os.path.join(input_dir, f) + ' -o ' + os.path.join(input_dir, f + '.pdf'))

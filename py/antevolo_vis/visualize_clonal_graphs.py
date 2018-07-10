#!/usr/bin/env python2

import os
import sys
import shutil

input_dir = sys.argv[1]
output_dir = sys.argv[2]
if os.path.exists(output_dir):
    shutil.rmtree(output_dir)
os.mkdir(output_dir)

files = os.listdir(input_dir)
for f in files:
    os.system('dot -Tpdf ' + os.path.join(input_dir, f) + ' -o ' + os.path.join(output_dir, f.split('.')[0] + '.pdf'))

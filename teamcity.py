#!/usr/bin/env python

import os.path
import sys
import shutil

def run_test():
    shutil.rmtree('build', True)
    exit_code = os.system('./prepare_cfg')
    if exit_code != 0:
        print('Preparing configuration files finished abnormally with exit code ' + str(exit_code))
        sys.exit(1)
    exit_code = os.system('make -j 8')
    if exit_code != 0:
        print('Compilation finished abnormally with exit code ' + str(exit_code))
        sys.exit(2)
    cmd = './ig_repertoire_constructor.py --test'
    exit_code = os.system(cmd)
    if exit_code != 0:
        print('IgRepCon finished abnormally with exit code ' + str(exit_code))
        sys.exit(3)
    output_dir = 'igrepcon_test'
    os.system('chmod -R 777 ' + output_dir)
    log = open(os.path.join(output_dir, 'ig_repertoire_constructor.log')).readlines()
    if len(log) < 2 or log[-2].strip() != 'Thank you for using IgRepertoireConstructor!':
        print 'Something goes wrong in IgRepertoireConstructor run'
        sys.exit(4)

if __name__ == '__main__':
    run_test()

#!/usr/bin/env python2

import os.path
import os
import argparse
import tempfile
import filecmp
import sys
import shutil


def random_id():
    import uuid
    return str(uuid.uuid4().get_hex().upper()[0:6])


def ParseCommandLine():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", "-i",
                        type=str,
                        help="Input reads for cropping",
                        required=True)
    parser.add_argument("--output", "-o",
                        type=str,
                        help="Desired output",
                        required=True)
    parser.add_argument("--args",
                        type=str,
                        default="",
                        help="Additional arguments for vj_finder")
    return parser.parse_args()


def InitParams(params):
    params.ig_kplus_vj_finder = os.path.join(os.path.dirname(os.path.realpath(__file__)),
        '../../../build/release/bin/ig_kplus_vj_finder')
    params.germline_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)),
        '../germline')


if __name__ == "__main__":
    print "Command line: %s" % " ".join(sys.argv)

    params = ParseCommandLine()
    InitParams(params)
    tmp_dir = tempfile.gettempdir()
    params.out_dir = "%s/%s" % (tmp_dir, random_id())

    cmd_line = "%(ig_kplus_vj_finder)s --db-directory=%(germline_dir)s --out=%(out_dir)s -i %(input)s %(args)s > /dev/null 2>&1" % params.__dict__
    os.system(cmd_line)

    if filecmp.cmp("%(out_dir)s/cropped.fa" % params.__dict__, params.output):
        print "Test passed"
    else:
        print "Test failed"

    shutil.rmtree(params.out_dir)

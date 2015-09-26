#!/usr/bin/env python2

import subprocess
import argparse


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Make test dataset using largest connectivity components")
    parser.add_argument("-n", "--count",
                        type=int,
                        default=13,
                        help="the number of joined connectivity components [default %(default)d]")
    parser.add_argument("-e", "--exclude-n-largest",
                        type=int,
                        default=0,
                        help="the number of excluded largest conn. comps. [default %(default)s]")
    parser.add_argument("input_dir",
                        type=str,
                        help="directory path with FASTQ files with connectivity components")
    parser.add_argument("output_file",
                        type=str,
                        help="output FASTQ file name")

    args = parser.parse_args()

    s = subprocess.check_output("wc -l %s/*.fastq | sort -n | tail -n %d | head -n %d" %
                                (args.input_dir, args.count + 1 + args.exclude_n_largest, args.count),
                                shell=True)

    l = list(s.split("\n"))
    l = [_.strip() for _ in l if _.strip()]
    l = [_.split(" ")[1] for _ in l]

    cmdline = "cat " + " ".join(l) + " > %s" % args.output_file
    subprocess.check_output(cmdline, shell=True)

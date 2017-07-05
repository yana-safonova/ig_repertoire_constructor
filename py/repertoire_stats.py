#!/usr/bin/env python2

from igquast_impl import get_clusters_sizes
from argparse import ArgumentParser

if __name__ == "__main__":
    parser = ArgumentParser(description="Report repertoire statistics")
    parser.add_argument("input",
                        type=str,
                        help="input repertoire file")

    args = parser.parse_args()

    print "Analized repertoire: %s" % args.input

    sizes = get_clusters_sizes(args.input)
    print "Reference consists of %d reads" % sum(sizes)
    print "Reference consists of %d clusters" % len(sizes)
    print "Reference consists of %d large (>=5) clusters" % len([size for size in sizes if size >= 5])
    print "Max cluster size %d" % max(sizes)

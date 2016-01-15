#!/usr/bin/env python2

import urllib
import os


def download_fasta(url, filename):
    f = urllib.urlopen(url)
    l = list(f)
    l = l[63:]
    l = [_.strip() for _ in l]
    i = 0
    while l[i]:
        i += 1
    l = l[:i]
    l = [_ + "\n" for _ in l]

    for i in xrange(len(l)):
        if l[i][0] != '>':
            l[i] = l[i].replace("y", "n").replace(".", "n")

    with open(filename, "w") as fh:
        fh.writelines(l)

eng2lat = {"human": "Homo sapiens",
           "rabbit": "Oryctolagus cuniculus",
           "mouse": "Mus",
           "rat": "Rattus norvegicus",
           "rainbow trout": "Oncorhynchus mykiss",
           "rhesus monkey": "Macaca mulatta",
           "pig": "Sus scrofa",
           "zebrafish": "Danio rerio",
           "platypus": "Ornithorhynchus anatinus",
           "bovine": "Bos taurus",
           "dog": "Canis lupus familiaris",
           "camel": "Camelus dromedarius"}
TCRs = ["TRAV", "TRAJ", "TRAC",
        "TRBV", "TRBD", "TRBJ", "TRBC",
        "TRGV", "TRGJ", "TRGC",
        "TRDV", "TRDD", "TRDJ", "TRDC"]

BCRs = ["IGHV", "IGHD", "IGHJ", "IGHC",
        "IGKV", "IGKJ", "IGKC",
        "IGLV", "IGLJ", "IGLC",
        "IGIV", "IGIJ", "IGIC"]


def download(locus, organism="human", allP=False, filename=""):
    assert locus in TCRs + BCRs
    assert organism in eng2lat.keys()

    organism_lat = eng2lat[organism]

    query_index = 2 if allP else 5

    url = "http://www.imgt.org/IMGT_GENE-DB/GENElect?query=7.%d+%s&species=%s" % (query_index, locus, organism_lat.replace(" ", "+"))
    if not filename:
        if allP:
            filename = "%s/%s-allP.fa" % (organism.replace(" ", "_"), locus)
        else:
            filename = "%s/%s.fa" % (organism.replace(" ", "_"), locus)

    download_fasta(url, filename)
    print "Organism: %s, locus: %s downloaded" % (organism, locus)

if __name__ == "__main__":
    for organism in ["human", "mouse"] + ["pig", "rat", "rabbit", "rhesus monkey"]:
        try:
            os.mkdir(organism.replace(" ", "_"))
        except:
            pass

    for locus in TCRs + BCRs:
        for organism in ["human", "mouse"]:
            download(locus, organism)
            download(locus, organism, True)

    for locus in BCRs:
        for organism in ["pig", "rat", "rabbit", "rhesus monkey"]:
            download(locus, organism)
            download(locus, organism, True)

import os
import sys

from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO

def ComputeNuclDict(align, pos):
    nucl_dict = dict()
    for i in range(0, len(align)):
        nucl = align[i][pos]
        if nucl == '-':
            continue
        if not nucl in nucl_dict:
            nucl_dict[nucl] = 0
        nucl_dict[nucl] += 1
    return nucl_dict

def NuclDictIsGood(nucl_dict, num_rows):
    num_cur_rows = 0
    for nucl in nucl_dict:
        num_cur_rows += nucl_dict[nucl]
    return float(num_cur_rows) / float(num_rows) > .5

def GetDominantPosition(nucl_dict):
    max_count = 0
    max_nucl = ""
    for nucl in nucl_dict:
        if nucl_dict[nucl] > max_count:
            max_count = nucl_dict[nucl]
            max_nucl = nucl
    return max_nucl, max_count

def GetRecessivePOsition(nucl_dict):
    dom_nucl, dom_count = GetDominantPosition(nucl_dict)
    max_count = 0
    max_nucl = ""
    for nucl in nucl_dict:
        if nucl_dict[nucl] > max_count and nucl != dom_nucl:
            max_count = nucl_dict[nucl]
            max_nucl = nucl
    return max_nucl, max_count

class AlignmentPosition:
    def __init__(self, dom_nucl, dom_count, rec_nucl, rec_count):
        self.dom_nucl = dom_nucl
        self.dom_count = dom_count
        self.rec_nucl = rec_nucl
        self.rec_count = rec_count

def main():
    if len(sys.argv) != 3:
        print "ERROR: invalid input arguments"
        print "python analyze_barcode_alignment.py barcode_alignment.aln output.txt"
        sys.exit(1)
    aln_fname = sys.argv[1]
    if not os.path.exists(aln_fname):
        print "ERROR: Alignment file " + aln_fname + " was not found"
        sys.exit(1)

    align = AlignIO.read(aln_fname, "clustal")
    num_columns = len(align[0])
    num_rows = len(align)
    print "Alignment for " + str(num_rows) + " rows & " + str(num_columns) + " columns was extracted from " + aln_fname
    amp_positions = list()
    for i in range(0, num_columns):
        nucl_dict = ComputeNuclDict(align, i)
        if not NuclDictIsGood(nucl_dict, num_rows):
            continue
        dom_nucl, dom_count = GetDominantPosition(nucl_dict)
        rec_nucl, rec_count = GetRecessivePOsition(nucl_dict)
        if float(rec_count) / float(dom_count) > .1:
            print str(i) + ": " + str(nucl_dict)
            print "Position " + str(i) + " looks like amplification"
            amp_positions.append(AlignmentPosition(dom_nucl, dom_count, rec_nucl, rec_count))
    print str(len(amp_positions)) + " amplification positions were computed from alignment file"

    output_fname = sys.argv[2]
    fhandler = open(output_fname, "w")
    for pos in amp_positions:
        fhandler.write(pos.dom_nucl + "\t" + str(pos.dom_count) + "\t" + pos.rec_nucl + "\t" + str(pos.rec_count) + "\n")
    print "Amplification positions were written to " + output_fname
    fhandler.close()

if __name__ == '__main__':
    main()
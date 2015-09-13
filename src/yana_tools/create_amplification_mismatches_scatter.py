import os
import sys

import numpy as np
import matplotlib.pyplot as plt

import operator
from matplotlib.ticker import FuncFormatter, MultipleLocator

class StatsRow:
    def __init__(self, file_row):
        splits = file_row.strip().split()
        if len(splits) != 4:
            print "ERROR: Invalid file row: " + file_row
            sys.exit(1)
        self.dom_nucl = splits[0]
        self.dom_count = splits[1]
        self.rec_nucl = splits[2]
        self.rec_count = splits[3]

def CreateMismatchesDict(file_rows):
    mismatches_dict = dict()
    for row in file_rows:
        mismatch_str = row.dom_nucl + row.rec_nucl
        if not mismatch_str in mismatches_dict:
            mismatches_dict[mismatch_str] = 0
        mismatches_dict[mismatch_str] += 1
    return mismatches_dict

def GetDigitByNucl(nucl):
    if nucl == 'A' or nucl == 'a':
        return 1
    if nucl == 'C' or nucl == 'c':
        return 2
    if nucl == 'G' or nucl == 'g':
        return 3
    if nucl == 'T' or nucl == 't':
        return 4
    print "ERROR: nucleotide " + nucl + " is unknown"
    sys.exit(1)

def GetNuclByDigit(digit, p = None):
    if digit == 1:
        return 'A'
    if digit == 2:
        return 'C'
    if digit == 3:
        return 'G'
    if digit == 4:
        return 'T'
    #print "ERROR: digit " + str(digit) + " is unknown"

def CreateXY(rec):
    nucl1 = rec[0]
    nucl2 = rec[1]
    return GetDigitByNucl(nucl1), GetDigitByNucl(nucl2)

def CreateXYArea(mismatches_dict):
    x = list()
    y = list()
    area = list()
    for rec in sorted(mismatches_dict):
        x_coord, y_coord = CreateXY(rec)
        x.append(x_coord)
        y.append(y_coord)
        area.append(mismatches_dict[rec] * np.pi)
    return x, y, area

def CreateScatterPlot(mismatches_dict):
    x, y, area = CreateXYArea(mismatches_dict)
    if len(x) != len(y) or len(x) != len(area):
        print "ERROR: invalid x, y, area lists"
        sys.exit(1)
    colors = np.random.rand(len(x))
    fig, ax = plt.subplots()
    plt.scatter(x, y, s = area, c = colors, alpha = 0.5)
    plt.xlabel("From ->")
    plt.ylabel("-> To")
    ax.xaxis.set_major_formatter(FuncFormatter(GetNuclByDigit))
    ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.yaxis.set_major_formatter(FuncFormatter(GetNuclByDigit))
    ax.yaxis.set_major_locator(MultipleLocator(1))
    plt.show()

def main():
    if len(sys.argv) != 2:
        print "ERROR: Invalid input parameters"
        print "python create_amplification_mismatches_scatter.py stats.txt"
        sys.exit(1)

    infile = sys.argv[1]
    if not os.path.exists(infile):
        print "ERROR: Input file " + infile + " was not found"
        sys.exit(1)

    ifhandler = open(infile, "r")
    rows = list()
    for l in ifhandler:
        l = l.strip()
        if l == "":
            break
        rows.append(StatsRow(l))
    print str(len(rows)) + " rows were extracted from " + infile

    mismatches_dict = CreateMismatchesDict(rows)
    print str(len(mismatches_dict)) + " mismatches pairs were extracted from " + infile
    for rec in sorted(mismatches_dict):
        print rec + ": " + str(mismatches_dict[rec])

    CreateScatterPlot(mismatches_dict)

if __name__ == '__main__':
    main()

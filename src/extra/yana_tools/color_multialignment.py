import os
import sys
from Bio import pairwise2

from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO

import operator

#colors = ['#5BC0EB', '#FDE74C', '#9BC53D', '#E55934', '#FA7921']
#colors = ["#EF597B", "#FF6D31", "#73B66B", "#FFCB18", "#29A2C6"]
colors = ["#5BC0EB", "#EF597B", "#FFCB18", "#FF6D31", "#73B66B"]


def CreateColumnDict(align, col_id):
    column_dict = dict()
    for i in range(0, len(align)):
        if align[i][col_id] not in column_dict:
            column_dict[align[i][col_id]] = 0
        column_dict[align[i][col_id]] += 1
    return column_dict

def AssignColorsDict(column_dict):
    color_dict = dict()
    sorted_column = sorted(column_dict.items(), key=operator.itemgetter(1), reverse=True)
    index = 0
    for ch in sorted_column:
        color_dict[ch[0]] = colors[index]
        index += 1
    return color_dict

def CreateSpaceDictBySeqName(align):
    max_rec_id = 0
    for rec in align:
        max_rec_id = max(max_rec_id, len(rec.id))
    name_n_space = max_rec_id + 4
    seq_dict = dict()
    for rec in align:
        seq_dict[rec.id] = name_n_space - len(rec.id)
    return seq_dict

def main(argv):
    if len(argv) != 3:
        print "Invalid input arguments"
        print "python color_multialignment.py alignment.align output.html"
        sys.exit(1)

    align_file = argv[1]
    if not os.path.exists(align_file):
        print "Alignment file " + align_file + " was not found"
        sys.exit(1)
    align = AlignIO.read(align_file, "clustal")
    num_columns = len(align[0])
    num_rows = len(align)
    print "Alignment (" + str(num_columns) + " columns, " + str(num_rows) + " rows) was extracted from " + align_file

    pos_color_dict = dict()
    for i in range(0, num_columns):
        column_dict = CreateColumnDict(align, i)
        if len(column_dict) == 1:
            continue
        color_dict = AssignColorsDict(column_dict)
        pos_color_dict[i] = color_dict
    print "Colors for positions were computed"

    html_fname = argv[2]
    ofhandler = open(html_fname, "w")
    ofhandler.write("<style>\n")
    ofhandler.write("    p {\n") 
    ofhandler.write("        font-family: monospace; \n")
    ofhandler.write("        font-size: 16px; \n")
    ofhandler.write("    }\n")
    ofhandler.write("</style>\n\n")

    ofhandler.write("<p>\n")    
    seq_dict = CreateSpaceDictBySeqName(align)
    chunk_size = 50
    num_blocks = num_columns / chunk_size
    if num_columns % chunk_size != 0:
        num_blocks += 1
    for i in range(0, num_blocks):
        start_pos = i * chunk_size
        end_pos = min(num_columns, (i + 1) * chunk_size)
        rec_ind = 0
        for rec in align:
            num_spaces = seq_dict[rec.id]
            ofhandler.write(rec.id + "&nbsp;" * num_spaces)
            for j in range(start_pos, end_pos):
                if j not in pos_color_dict:
                    ofhandler.write(align[rec_ind][j])
                else:
                    color = pos_color_dict[j][align[rec_ind][j]]
                    ofhandler.write("<font style=\"BACKGROUND-COLOR:" + color + "\">" + align[rec_ind][j] + "</font>")
            ofhandler.write("<br>")
            rec_ind += 1
        ofhandler.write("<br><br>")
    ofhandler.write("</p>")
    ofhandler.close()
    print "Colored alignmnent was written to " + html_fname

if __name__ == '__main__':
    main(sys.argv)


import os
import sys
sys.path.append('py/utils/')
import drawing_utils

if len(sys.argv) != 2:
    print "nvalid params"
    sys.exit(1)
fname = sys.argv[1] #'opt_kmers_tau_4_k_35.txt'
lines = open(fname, 'r').readlines()
rel_pos = list()
for l in lines:
    splits = l.strip().split()
    length = int(splits[1])
    for i in range(2, len(splits)):
        rel_pos.append(float(splits[i]) / length)
#print rel_pos
figure_fname = fname.split('.')[0] + '.pdf'
settings = drawing_utils.GetGraphicalSettings(xlabel = "Representative k-mer positions", ylabel = "# k-mers", title = "", output_filename = figure_fname, bins = 100, label = "", histtype = "bar", xlog_scale = False, ylog_scale = False, draw_legend = False, colors = "", legend_loc = 'upper right', show_xaxis = True, show_yaxis = True, xmin_shift = 0, xmax_shift = 0, ymin_shift = 0, ymax_shift = 0)
drawing_utils.DrawMutationHistogram(rel_pos, settings)

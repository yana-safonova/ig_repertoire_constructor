import logging
import sys

import matplotlib as mplt
mplt.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def get_logger_by_arg(logger_arg, logger_name = ""):
    if logger_arg != "-":
        return logger_arg
    log = logging.getLogger(logger_name)
    log.setLevel(logging.DEBUG)
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter('%(message)s'))
    console.setLevel(logging.DEBUG)
    log.addHandler(console)
    return log

def output_figure(output_fname, log_info, log):
    pp = PdfPages(output_fname + ".pdf")
    pp.savefig()
    pp.close()
    plt.savefig(output_fname + ".png")
    plt.savefig(output_fname + ".svg")
    plt.clf()
    log.info(log_info + " was written to " + output_fname + ".pdf, .png and .svg")

hydrophobic = ['L', 'I', 'F', 'W', 'V', 'M', 'C', 'Y', 'A']
neutral = ['T', 'E', 'G', 'S', 'Q', 'D', '*']
hydrophilic = ['R', 'K', 'N', 'H', 'P']

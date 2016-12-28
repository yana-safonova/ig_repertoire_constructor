#!/usr/bin/env python2

from simulate import *


from time import time



if __name__ == "__main__":
    tmp_dir = igrec_dir + "/timing_tmp_dir"
    igrec_times = []
    mixcr_times = []
    vjf_times = []
    kaligner_times = []
    for i in xrange(6):
    # for i in xrange(1):
        filename = igrec_dir + "SRA_performance_test/cleaned_reads/SRR138346%d.fa.gz" % i
        # filename = igrec_dir + "SRA_performance_test/merged_reads/SRR138346%d.fq.gz" % i

        print "testing " + filename
        real_time = time()
        run_igrec(filename, tmp_dir, tau=3, additional_args=" --no-alignment")
        real_time = time() - real_time
        igrec_times.append(real_time)

        real_time = time()
        run_mixcr2(filename, tmp_dir)
        real_time = time() - real_time
        mixcr_times.append(real_time)

        real_time = time()
        run_vjfinder(filename, tmp_dir)
        real_time = time() - real_time
        vjf_times.append(real_time)

        real_time = time()
        run_mixcr2_alignment_only(filename, tmp_dir)
        real_time = time() - real_time
        kaligner_times.append(real_time)

    for times in [igrec_times, mixcr_times, vjf_times, kaligner_times]:
        print sum(times), times

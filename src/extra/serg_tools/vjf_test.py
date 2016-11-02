from Bio import SeqIO


def main():
    reads_orig = []
    reads_orig.extend([read.name for read in SeqIO.parse("big_cluster.fastq", "fastq")])
    reads_orig.extend([read.name for read in SeqIO.parse("added_reads.fastq", "fastq")])
    print len(reads_orig)

    reads = [read for read in SeqIO.parse("/home/serg/local_data/igrec_umi/vj_finder/cleaned_reads.fa", "fasta") if read.name in reads_orig]
    print len(reads)
    for read1 in reads:
        s1 = str(read1.seq)
        for read2 in reads:
            s2 = str(read2.seq)
            if s1 != s2 and s1 in s2:
                print read1.name
                print read1.seq
                print read2.name
                print read2.seq
                exit(0)
                print "-----------"


if __name__ == '__main__':
    main()

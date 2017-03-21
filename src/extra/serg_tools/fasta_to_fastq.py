def main():
    from Bio import SeqIO
    import sys
    import os
    # with open(sys.argv[1]) as fasta_file, open(os.path.splitext(sys.argv[1])[0] + '.fastq', "w") as fastq_file:
        # records = SeqIO.parse(fasta_file, "fasta")
        # for record in records:
        #     record.qual = "H" * len(record.seq)
        # SeqIO.write(records, fastq_file, "fastq")


if __name__ == '__main__':
    main()

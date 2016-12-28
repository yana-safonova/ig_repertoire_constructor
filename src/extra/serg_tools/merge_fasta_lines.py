def main():
    import sys
    import os
    fname_ext = os.path.splitext(sys.argv[1])
    first = True
    with open(sys.argv[1]) as fasta_file, open(fname_ext[0] + '_1' + fname_ext[1], "w") as output_file:
        for line in fasta_file:
            if line[0] == '>':
                if not first:
                    output_file.write("\n")
                output_file.write(line)
            else:
                output_file.write(line.strip())
            first = False


if __name__ == '__main__':
    main()

def ParseArgs():
    from argparse import ArgumentParser
    parser = ArgumentParser(description="Compress equal clusters")
    parser.add_argument("input",
                        type=str,
                        help="input FASTA file")
    parser.add_argument("output",
                        type=str,
                        help="output FASTA file")
    parser.set_defaults(count_barcodes=False)

    return parser.parse_args()


def main():
    args = ParseArgs()
    with open(args.input) as in_file:
        with open(args.output, "w") as out_file:
            for line in in_file:
                if line[0] == '>':
                    import re
                    match = re.match(r"^>antibody_(\d+)_multiplicity_(\d+)_copy_\d+$", line)
                    groups = match.groups()
                    out_file.write(">cluster___%d___size___1\n" % int(groups[0]))
                else:
                    out_file.write(line)


if __name__ == "__main__":
    main()

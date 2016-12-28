import os
from Bio import SeqIO
import random
import numpy as np
from Bio import Seq


# >antibody_1_multiplicity_3_copy_1
def parse_header(s):
    s = s.strip().split("_")
    assert s[0] == "antibody"
    assert s[2] == "multiplicity"
    assert s[4] == "copy"

    return map(int, s[1::2])


if __name__ == "__main__":
    res = 0

    with open("final_repertoire.fa") as f:
        for record in SeqIO.parse(f, "fasta"):
            name = str(record.description)
            id, mult, copy = parse_header(name)
            if copy == 1 and mult >= 5:
                res += 1

    print res

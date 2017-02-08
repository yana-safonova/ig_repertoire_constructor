import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input",
                        help="Config with input files",
                        required=True)
    return parser.parse_args()

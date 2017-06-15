import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input",
                        help="Config with input files",
                        # required=True,
                        default="/Sid/abzikadze/datasets/config.json")
    parser.add_argument("--skip-analysis",
                        action="store_const",
                        const=True,
                        dest="skip_analysis",
                        help="Only estimate model")
    return parser.parse_args()

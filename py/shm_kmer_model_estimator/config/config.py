import json
import os.path
from argparse import Namespace

current_dir = os.path.dirname(os.path.realpath(__file__))

configfile = os.path.join(current_dir, "config.json")
with open(configfile, "r") as handle:
    config = json.load(handle, object_hook=lambda d: Namespace(**d))

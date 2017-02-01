import json
import os.path
from argparse import Namespace


configfile = os.path.join("config", "config.json")
with open(configfile, "r") as handle:
    config = json.load(handle, object_hook=lambda d: Namespace(**d))

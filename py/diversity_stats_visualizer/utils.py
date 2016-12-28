import logging
import sys

def get_logger_by_arg(logger_arg, logger_name = ""):
    if logger_arg != "-":
        return logger_arg
    log = logging.getLogger(logger_name)
    log.setLevel(logging.DEBUG)
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter('%(message)s'))
    console.setLevel(logging.DEBUG)
    log.addHandler(console)
    return log
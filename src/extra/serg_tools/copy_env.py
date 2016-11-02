import fnmatch
import os
import shutil
import sys


def make_dir_if_not_there(dest):
    try:
        os.makedirs(dest)
    except OSError:
        if not os.path.isdir(dest):
            raise


def main():
    dest = sys.argv[1]
    clean = int(sys.argv[2]) > 0 if len(sys.argv) > 2 else True
    if clean:
        shutil.rmtree(dest, ignore_errors = True)
    make_dir_if_not_there(dest)
    current_dir = os.path.dirname(os.path.realpath(__file__))
    igrec_dir = os.path.join(current_dir, os.pardir, os.pardir, os.pardir)
    for root, dirs, files in os.walk(igrec_dir):
        for file in fnmatch.filter(files, "*.py"):
            dir = os.path.relpath(root, igrec_dir)
            dest_dir = os.path.join(dest, dir)
            make_dir_if_not_there(dest_dir)
            shutil.copy2(os.path.join(root, file), dest_dir)
    for dir in ["build", "pipeline_makefiles", "configs", "data"]:
        shutil.copytree(os.path.join(igrec_dir, dir), os.path.join(dest, dir))

if __name__ == '__main__':
    main()

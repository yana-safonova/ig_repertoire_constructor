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


def copy_file(dest, file, igrec_dir, root):
    dir = os.path.relpath(root, igrec_dir)
    dest_dir = os.path.join(dest, dir)
    make_dir_if_not_there(dest_dir)
    # print "copying %s to %s" % (os.path.join(root, file), dest_dir)
    shutil.copy2(os.path.join(root, file), dest_dir)


def exclude_broken_links(dir, files):
    return [file for file in files if os.path.splitext(file)[1] == 'py']


def main():
    dest = sys.argv[1]
    clean = int(sys.argv[2]) > 0 if len(sys.argv) > 2 else True
    if clean:
        shutil.rmtree(dest, ignore_errors = True)
    make_dir_if_not_there(dest)
    current_dir = os.path.dirname(os.path.realpath(__file__))
    igrec_dir = os.path.join(current_dir, os.pardir, os.pardir, os.pardir)
    walk = list(os.walk(igrec_dir))
    for root, dirs, files in walk:
        for file in [file for file in files if file.endswith((".py", ".jar", ".sh")) and '/build/' not in root and '/cmake-' not in root]:
            print "copying %s/%s" % (root, file)
            copy_file(dest, file, igrec_dir, root)
    for file in ["test_dataset/merged_reads.fastq"]:
        copy_file(dest, os.path.split(file)[-1], igrec_dir, os.path.join(igrec_dir, os.path.dirname(file)))
    for dir in ["build", "pipeline_makefiles", "configs", "data"]:
        shutil.copytree(os.path.join(igrec_dir, dir), os.path.join(dest, dir), exclude_broken_links)

    os.system("git -C %s log | head > %s/GIT_REVISION" % (igrec_dir, dest))


if __name__ == '__main__':
    main()

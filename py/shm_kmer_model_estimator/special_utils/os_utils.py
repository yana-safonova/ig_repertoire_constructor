import os
import errno

def list_only_dirs(abs_path):
    return filter(lambda x: os.path.isdir(os.path.join(abs_path, x)),
                            os.listdir(abs_path))


def smart_mkdir(dirname):
    try:
        os.mkdir(dirname)
    except OSError as exc:
        if exc.errno != errno.EEXIST:
            raise exc


def smart_makedirs(dirname):
    try:
        os.makedirs(dirname)
    except OSError as exc:
        if exc.errno != errno.EEXIST:
            raise exc


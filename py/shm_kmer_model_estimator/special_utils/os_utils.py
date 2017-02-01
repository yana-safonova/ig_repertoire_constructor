import os

def list_only_dirs(abs_path):
    return filter(lambda x: os.path.isdir(os.path.join(abs_path, x)),
                            os.listdir(abs_path))

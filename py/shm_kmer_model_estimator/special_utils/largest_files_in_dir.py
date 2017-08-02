""" author: http://codereview.stackexchange.com/questions/8958/n-largest-files-in-a-directory """
import os
import heapq


def n_largest_files(directory, num_files=5):
    file_names = \
        (os.path.join(path, name) for path, _, filenames in os.walk(directory)
         for name in filenames)
    big_files = heapq.nlargest(num_files, file_names, key=os.path.getsize)
    return big_files

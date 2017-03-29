#!/usr/bin/env python

import os
import sys
import subprocess

class ThreadException(Exception):
    code = 0
    def __init__(self, value, _code):
        self.value = value
        self.code = _code

class JobQueue:
    worker_list = []
    def __init__(self, threads):
        self.worker_list = [None] * threads

    def interrupt(self):
        for worker in self.worker_list:
            if worker != None:
                try:
                    worker.kill()
                except OSError:
                    pass


    def find_worker(self):
        while True:
            for i in range(len(self.worker_list)):
                worker = self.worker_list[i]
                if worker == None:
                    return i
                if not worker.poll() and worker.returncode is None:
                    continue
                if worker.returncode is not None and worker.returncode != 0:
                    raise ThreadException('Thread finished abnormally', worker.returncode)
                worker = None
                return i

    def wait(self):
        for worker in self.worker_list:
            if worker is not None:
                worker.wait()
                if worker.returncode != 0:
                    raise ThreadException('Thread finished abnormally', worker.returncode)

    def run(self, job_list):
        try:
            for job in job_list:
                worker_id = self.find_worker()
                self.worker_list[worker_id] = subprocess.Popen(job, shell = True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, stdin=subprocess.PIPE)
            self.wait()
        except ThreadException as e:
            self.interrupt()
            return e.code
        return 0

def SimulateXargs(parameter_list, command_line, threads):
    command_list = [command_line.replace("{}", parameter) for parameter in parameter_list]
    return JobQueue(threads).run(command_list)

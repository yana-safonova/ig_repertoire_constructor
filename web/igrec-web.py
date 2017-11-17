#!/usr/bin/env python2

# TODO
# nice file uploader
# https://www.startutorial.com/articles/view/resumable-file-upload-part-1
# https://github.com/23/resumable.js/blob/master/samples/Frontend%20in%20jQuery.md


from flask import Flask
from flask import render_template
from flask import request
from flask import make_response
from flask import redirect
from flask import url_for
from flask import send_from_directory
from flask_autoindex import AutoIndex
from flask import session
from flask_session import Session

from flask_socketio import SocketIO


app = Flask(__name__)


import os
current_dir = os.path.dirname(os.path.realpath(__file__))
igrec_dir = os.path.dirname(current_dir)

runs_dir = current_dir + "/runs"
idx = AutoIndex(app, runs_dir, add_url_rules=False)


app.config['SESSION_TYPE'] = "filesystem"  # TODO Use redis here as well
app.config['SESSION_USE_SIGNER'] = True
app.config['SESSION_FILE_DIR'] = os.path.join(current_dir, "sessions")
app.config['SECRET_KEY'] = "dc6627ce-2b94-4c45-a1e7-75fd595e2aca"

app.config['CELERY_RESULT_BACKEND'] = "redis://localhost:6379/0"
app.config['CELERY_BROKER_URL'] = "redis://localhost:6379/0"  # FIXME Use user-defined db-id
# app.config['CELERY_BROKER_URL'] = "amqp://user:password@localhost:5672/vhost"  # FIXME Use user-defined db-id
# TODO Set CELERY_TRACK_STARTED


# app.config["USE_X_SENDFILE"] = True # TODO Add nginx support

import redis
r = redis.StrictRedis(host='localhost', port=6379, db=0)

Session(app)

socketio = SocketIO(app)

from celery import Celery

def make_celery(app):
    celery = Celery(app.import_name, backend=app.config['CELERY_RESULT_BACKEND'],
                    broker=app.config['CELERY_BROKER_URL'])
    celery.conf.update(app.config)
    TaskBase = celery.Task
    class ContextTask(TaskBase):
        abstract = True
        def __call__(self, *args, **kwargs):
            with app.app_context():
                return TaskBase.__call__(self, *args, **kwargs)
    celery.Task = ContextTask
    return celery


celery = celery_app = make_celery(app)

@celery_app.task()
def run_command_simple(cmd):
    import os
    rc = os.system(cmd)
    return rc



def terminate_task(id, signal="SIGTERM"):
    from celery import Terminated
    try:
        celery.control.revoke(id, terminate=True, signal=signal)
    except Terminated:
        pass


@app.route('/tasks/commands/terminate/<uuid:task_id>')
def tasks_command_terminate(task_id):
    terminate_task(str(task_id))
    return "Terminated"


@app.route('/tasks/commands/get_status/<uuid:task_id>')
def tasks_command_get_status(task_id):
    task = execute.AsyncResult(str(task_id))

    return task.status

@app.route("/tasks/commands/get_output/<uuid:task_id>")
def tasks_commands_get_output(task_id):
    task = execute.AsyncResult(str(task_id))

    if task.status in ["SUCCESS"]:
        return task.info["output_id"]
    else:
        return ""


@celery_app.task(bind=True)
def execute(self, command, output_id):
    import subprocess
    from celery.exceptions import Ignore

    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    log = []
    while process.poll() is None:
        nextline = process.stdout.readline()
        log.append(nextline)
        self.update_state(state="PROGRESS", meta={"log": log})

    stdout, stderr = process.communicate()
    rc = process.returncode

    self.update_state(state="SUCCESS", meta={"log": log, "rc": rc,
                                             "stdout": stdout, "stderr": stderr,
                                             "output_id": output_id})

    # Setup state and meta manually
    # from https://stackoverflow.com/questions/25826639/how-to-manually-mark-a-celery-task-as-done-and-set-its-result
    raise Ignore()


def add_task(command, output_id):
    task = execute.delay(command, output_id)

    r.lpush("TaskList", task.id)
    return task






# TODO Use redis to real-time tracking
# https://stackoverflow.com/questions/35437668/realtime-progress-tracking-of-celery-tasks

@app.route("/tasks/status/<uuid:task_id>")
def status_page(task_id):
    task = execute.AsyncResult(str(task_id))

    if task.status in ["PROGRESS", "SUCCESS"]:
        log = task.info["log"]
    else:
        log = None
    if task.status in ["SUCCESS"]:
        output_id=task.info["output_id"]
    else:
        output_id=None

    return render_template("status.html", status=task.status, log=log,
                           task_id=task_id,
                           output_id=output_id)




@app.route("/jQueryFileTree", methods=["POST"])
def jQueryFileTree():
    from os import listdir
    from os.path import isfile, isdir, join

    path = request.form["dir"]
    try:
        files = [f for f in listdir(path) if isfile(join(path, f))]
        dirs = [f for f in listdir(path) if isdir(join(path, f))]
    except:
        # No access
        files = dirs = []

    # Filter hidden files
    def filter_out_hidden_files(filenames):
        return [f for f in filenames if not f.startswith(".")]

    files = filter_out_hidden_files(files)
    dirs = filter_out_hidden_files(dirs)

    def filter_fastx_files(filenames):
        import re
        return [f for f in files if re.match(r"^.*\.f(ast)?(a|q)(\.gz|\.bz2)?$", f)]

    files = filter_fastx_files(files)

    # TODO Check sorting options of jQueryFileTree module
    return render_template("jQueryFileTree.html", path=path,
                           files=sorted(files), dirs=sorted(dirs))


@app.route('/')
def new_run():
    response = make_response(render_template("index.html"))
    # response.set_cookie('username', 'the username')
    return response


@app.route('/report/<output_id>')
def report(output_id):
    response = make_response(render_template("report.html", output_id=output_id))
    # response.set_cookie('username', 'the username')
    return response


def create_uuid_dir(basedir):
    import uuid
    import os

    name = str(uuid.uuid4())
    try:
        os.mkdir(os.path.join(basedir, name))
    except OSError:
        return create_uuid_dir(basedir)

    return name


@app.route('/runs/')
@app.route('/runs/<path:path>')
def autoindex(path='.'):
    return idx.render_autoindex(path, browse_root=runs_dir)

# TODO
# Implement file download resuming as it shown there: http://flask.pocoo.org/docs/0.10/api/#flask.Flask.use_x_sendfile
# https://github.com/pallets/flask/issues/1388



@app.route('/run_igrec', methods=["POST"])
def run_igrec():
    form = dict(request.form)
    # TODO Add validation
    for k in form:
        form[k] = form[k][0]


    form["output_id"] = output_id = create_uuid_dir(runs_dir)
    form["output"] = output = os.path.join(runs_dir, output_id)
    form["igrec_dir"] = igrec_dir

    form["input_line"] = "-s %(merged-reads-file)s" % form if form["readstype"] == "merged" else "-1 %(paired-reads-file-right)s -2 %(paired-reads-file-left)s" % form

    form["tool"] = "igrec.py" if form["barcodingtype"] == "no" else "barcoded_igrec.py"
    if form["barcodingtype"] == "file":
        return "Barcodes in a separate file are not supported yet =("

    cmd = "%(igrec_dir)s/%(tool)s %(input_line)s \
        --loci=%(loci)s --organism=%(organism)s \
        --threads=3 --output=%(output)s" % form

    # print form
    if "skip-alignment" in form:
        cmd += " --no-alignment"

    if "do-igquast" in form:
        cmd += " && %(igrec_dir)s/igquast.py \
            --initial-reads=%(output)s/vjfinder/cleaned_reads.fa \
            --constructed-repertoire=%(output)s/final_repertoire.fa \
            --constructed-rcm=%(output)s/final_repertoire.rcm \
            --output=%(output)s/igquast \
            --reference-free" % form

    # print cmd
    cmd = "ping ya.ru -c 1000; " + cmd
    task = add_task(cmd, output_id=output_id)

    # task = execute.delay("ping ya.ru -c 15", output_id=output_id)

    session['last_output'] = output
    session['last_task'] = task

    return redirect(url_for("status_page", task_id=task.id))


@app.route('/tasks')
def task_list():
    tasks = r.lrange("TaskList", 0, -1)

    return render_template("jobs.html", tasks=tasks)


if __name__ == "__main__":
    socketio.run(app)

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
    from celery.exceptions import Terminated
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
    return get_task_status(task_id)

@app.route("/tasks/commands/get_output/<uuid:task_id>")
def tasks_commands_get_output(task_id):
    task_id = str(task_id)
    return r.get("output." + task_id) or ""


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
    import datetime

    now = datetime.datetime.now()
    task = execute.delay(command, output_id)

    r.rpush("TaskList", task.id)
    r.set("start_time." + task.id, now.strftime("%Y-%m-%d %H:%M:%S"))
    r.set("command." + task.id, command)
    r.set("output." + task.id, output_id)
    return task


@app.route('/admin')
def admin():
    return render_template("admin.html")



def rerun(task_id):
    command = r.get("command." + task_id)
    output = r.get("output." + task_id)
    if (command and output):
        new_output = create_uuid_dir(runs_dir, "_(rerun_%s)" % output)
        new_task = add_task(command, new_output)
        return new_task.id
    else:
        return False


@app.route('/tasks/commands/rerun/<uuid:task_id>')
def tasks_command_rerun(task_id):
    task_id = str(task_id)
    result = rerun(task_id)
    return result or "FAILED"



def get_task_status(task_id):
    task = execute.AsyncResult(str(task_id))

    status = task.status
    if task.status in ["SUCCESS"]:
        if task.info["rc"] != 0:
            status = "CRASHED"

    return status

# TODO Use redis to real-time tracking
# https://stackoverflow.com/questions/35437668/realtime-progress-tracking-of-celery-tasks



@app.route("/tasks/status/<uuid:task_id>")
def status_page(task_id):
    task_id = str(task_id)
    task = execute.AsyncResult(task_id)

    status = get_task_status(task_id)

    if status in ["PROGRESS", "SUCCESS", "CRASHED"]:
        log = task.info["log"]
    else:
        log = None

    output_id=r.get("output." + task_id)

    return render_template("status.html", status=status, log=log,
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


def create_uuid_dir(basedir, suffix=""):
    import uuid
    import os
    import time
    import hashlib

    id = str(uuid.uuid4())

    name = time.strftime("%d%b%Y_%H-%M-%S")
    name += "-" + hashlib.md5(id).hexdigest()[:5]
    name += suffix
    try:
        os.mkdir(os.path.join(basedir, name))
    except OSError:
        return create_uuid_dir(basedir, suffix=suffix)

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
        --threads=%(threads)s --output=%(output)s" % form
    if form["barcodingtype"] == "no":
        cmd += " --tau=%(tau)s" % form

    # print form
    if "skip-alignment" in form:
        cmd += " --no-alignment"


    if "skip-alignment" in form:
        if form["readstype"] == "merged":
            form["initial-reads-file"] = form["merged-reads-file"]
        else:
            form["initial-reads-file"] = "%(output)s/merged_reads.fq" % form
    else:
        form["initial-reads-file"] = "%(output)s/vj_finder/cleaned_reads.fa" % form

    if "do-igquast" in form:
        if form["barcodingtype"] == "no":
            cmd += " && %(igrec_dir)s/igquast.py \
                --initial-reads=%(initial-reads-file)s \
                --constructed-repertoire=%(output)s/final_repertoire.fa \
                --constructed-rcm=%(output)s/final_repertoire.rcm \
                --output=%(output)s/igquast \
                --reference-free" % form
        else:
            cmd += " && %(igrec_dir)s/igquast.py \
                --initial-reads=%(initial-reads-file)s \
                --constructed-repertoire=%(output)s/final_repertoire/final_repertoire.fa \
                --constructed-rcm=%(output)s/final_repertoire/final_repertoire.rcm \
                --output=%(output)s/igquast \
                --reference-free" % form

    print cmd
    # cmd = "ping ya.ru -c 1000; " + cmd
    task = add_task(cmd, output_id=output_id)

    # task = execute.delay("ping ya.ru -c 15", output_id=output_id)

    session['last_output'] = output
    session['last_task'] = task

    return redirect(url_for("status_page", task_id=task.id))


@app.route('/tasks/commands/purge')
def reset_all_tasts():
    # TODO Add unique command token
    for key in r.scan_iter("*"):
        print key
        r.delete(key)

    return "OK"


@app.route('/tasks')
def task_list():
    tasks = r.lrange("TaskList", 0, -1)

    def get_task_time(task_id):
       time = r.get("start_time." + task_id)
       return time

    times = map(get_task_time, tasks)

    return render_template("tasks.html", tasks=tasks, times=times, zipped=zip(tasks, times))

@app.route('/uploads')
def uploads():
    return render_template("upload.html")


# From here https://github.com/szelcsanyi/resumable.js/blob/8872d84c57b8c8cc756f7b74ec51b78233f82021/samples/Backend%20on%20Flask.md

temp_base = current_dir + "/tmp/"
upload_dir = current_dir + "/uploads/"
from flask import abort

@app.route('/upload', methods=['GET'])
def show():
    identifier = request.args.get('resumableIdentifier', type=str)
    filename = request.args.get('resumableFilename', type=str)
    chunk_number = request.args.get('resumableChunkNumber', type=int)

    if not identifier or not filename or not chunk_number:
        # Parameters are missing or invalid
        abort(500, 'Parameter error')

    # chunk folder path based on the parameters
    temp_dir = "{}{}".format(temp_base, identifier)
    # chunk path based on the parameters
    chunk_file = "{}/{}.part{}".format(temp_dir, filename, chunk_number)

    app.logger.debug('Getting chunk: %s', chunk_file)

    # TODO Get rid of code duplication
    if os.path.isfile(chunk_file):
        # Let resumable.js know this chunk already exists
        return 'OK'
    else:
        # Let resumable.js know this chunk does not exists and needs to be uploaded
        abort(404, 'Not found')



@app.route('/upload', methods=['POST'])
def create():
    identifier = request.form.get('resumableIdentifier', type=str)
    filename = request.form.get('resumableFilename', type=str)
    chunk_number = request.form.get('resumableChunkNumber', type=int)
    chunk_size = request.form.get('resumableChunkSize', type=int)
    current_chunk_size = request.form.get('resumableCurrentChunkSize', type=int)
    total_size = request.form.get('resumableTotalSize', type=int)

    if not identifier or not filename or not chunk_number or not chunk_size or not current_chunk_size or not total_size:
        # Parameters are missing or invalid
        abort(500, 'Parameter error')

    # chunk folder path based on the parameters
    temp_dir = "{}{}".format(temp_base, identifier)
    # chunk path based on the parameters
    chunk_file = "{}/{}.part{}".format(temp_dir, filename, chunk_number)

    app.logger.debug('Creating chunk: %s', chunk_file)

    # Create directory of not exists
    try:
        os.makedirs(temp_dir, 0777)
    except:
        pass

    file = request.files['file']
    file.save(chunk_file)

    from os import listdir
    from os.path import isfile, join
    import math
    onlyfiles = [f for f in listdir(temp_dir) if isfile(join(temp_dir, f))]

    # When all chunks are uploaded
    n_of_chunks = int(math.ceil(float(total_size) / chunk_size))
    print "#chunks", n_of_chunks
    if len(onlyfiles) == n_of_chunks:
        # Create a target file
        target_file_name = "{}/{}".format(upload_dir, filename)
        with open(target_file_name, "wb") as target_file:  # TODO Use flock here
            # Loop through the chunks
            for i in range(1, n_of_chunks + 1):
                # Select the chunk
                stored_chunk_file_name = "{}/{}.part{}".format(temp_dir, filename, str(i))
                with open(stored_chunk_file_name, 'rb') as stored_chunk_file:
                    # Write chunk into target file
                    target_file.write(stored_chunk_file.read())  # We do not need buffering here, chunks are already small enough
                # Deleting chunk
                os.unlink(stored_chunk_file_name)
        os.rmdir(temp_dir)
        app.logger.debug('File saved to: %s', target_file_name)

    return 'OK'




if __name__ == "__main__":
    socketio.run(app)

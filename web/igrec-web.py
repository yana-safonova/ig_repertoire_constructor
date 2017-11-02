#!/usr/bin/env python2

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


app.config['SESSION_TYPE'] = "filesystem"
app.config['SESSION_USE_SIGNER'] = True
app.config['SESSION_FILE_DIR'] = os.path.join(current_dir, "sessions")
app.config['SECRET_KEY'] = "dc6627ce-2b94-4c45-a1e7-75fd595e2aca"
Session(app)

socketio = SocketIO(app)


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
def index():
    response = make_response(render_template("index.html"))
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
    input = request.form["merged_reads"]
    output_id = create_uuid_dir(runs_dir)
    output = os.path.join(runs_dir, output_id)
    return_code = os.system("%s/igrec.py -s %s --loci=all --threads=3 --output=%s" % (igrec_dir, input, output))

    session['last_input'] = input
    session['last_output'] = output

    print return_code
    return redirect(url_for("autoindex", path=output_id))


if __name__ == "__main__":
    socketio.run(app)

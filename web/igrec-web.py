from flask import Flask, render_template
from flask import request
from flask import make_response
from flask import redirect
from flask import url_for
from flask import send_from_directory
from flask_autoindex import AutoIndex


import os
current_dir = os.path.dirname(os.path.realpath(__file__))
igrec_dir = os.path.dirname(current_dir)

app = Flask(__name__)
runs_dir = current_dir + "/runs"
idx = AutoIndex(app, runs_dir, add_url_rules=False)


def install_secret_key(app, filename='secret_key'):
    import sys
    import os.path
    """Configure the SECRET_KEY from a file
    in the instance directory.

    If the file does not exist, print instructions
    to create it from a shell with a random key,
    then exit.

    Taken from http://flask.pocoo.org/snippets/104/

    """
    filename = os.path.join(app.instance_path, filename)
    try:
        app.config['SECRET_KEY'] = open(filename, 'rb').read()
    except IOError:
        print 'Error: No secret key. Create it with:'
        if not os.path.isdir(os.path.dirname(filename)):
            print 'mkdir -p', os.path.dirname(filename)
        print 'head -c 24 /dev/urandom >', filename
        sys.exit(1)

install_secret_key(app)


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

    # TODO Check sorting options of jQueryFileTree module
    return render_template("jQueryFileTree.html", path=path,
                           files=sorted(files), dirs=sorted(dirs))


@app.route('/')
def index():
    response = make_response(render_template("index.html"))
    # response.set_cookie('username', 'the username')
    return response

# @app.route('/result/<output_id>')
# def result(output_id):
#     with open(current_dir + "/" + output_id + "/igrec.log") as log:
#         return render_template("log.html", lines=log)



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


@app.route('/run_igrec', methods=["POST"])
def run_igrec():
    input = request.form["merged_reads"]
    output_id = create_uuid_dir(runs_dir)
    output = os.path.join(runs_dir, output_id)
    return_code = os.system("%s/igrec.py -s %s --loci=all --threads=3 --output=%s" % (igrec_dir, input, output))

    print return_code
    return redirect(url_for("autoindex", path=output_id))



@app.route('/hello/')
@app.route('/hello/<name>')
def hello(name=None):
    response = make_response(render_template('hello.html', name=name))
    response.set_cookie('username', 'the username')
    return response

if __name__ == "__main__":
    app.run()

#!/usr/bin/env python2

import os
import os.path
import fnmatch
import sys
import json

current_dir = os.path.dirname(os.path.realpath(__file__))
igrec_dir = current_dir + "/../../../../"
sys.path.append(igrec_dir + "/py/")
from ash_python_utils import mkdir_p

def metadatajson():
    json_file = json.dumps(
        {
            "Name": "",
            "Description": "",
            "HrefAppSession": "",
            "Properties": [
                {
                    "Type": "sample[]",
                    "Name": "Input.Samples",
                    "Items": [

                    ]
                }
            ]
        }
    )
    return json_file


def showDir(dir):
    import os
    import fnmatch

    for root, dir, files in os.walk(dir):
        print root
        print ""
        for items in fnmatch.filter(files, "*"):
            print "..." + items
            print ""


def safe_cast(val, to_type, default=None):
    try:
        return to_type(val)
    except ValueError:
        return default


if __name__ == "__main__":
    jsonfile = open('/data/input/AppSession.json')
    jsonObject = json.load(jsonfile)

    runID = jsonObject['Id']
    runName = jsonObject['Name']
    runName = runName.replace(' ', '_')

    param = {}
    param["threads"] = 32
    param["min-cluster-size"] = 5
    param["tau"] = 4
    param["min-fillin"] = 0.6
    param["loci"] = "all"
    param["organism"] = "human"
    param["additionalfiles"] = False
    param["pseudogenes"] = False
    param["AppSessionName"] = "IgReC"

    # loop over properties
    for item in jsonObject['Properties']['Items']:
        name = item['Name']

        if name == 'Input.reads-file':
            ParentAppResult = item['Content']['ParentAppResult']
            sampleDir = "/data/input/appresults/%s" % ParentAppResult['Id']

            param['reads'] = "%s/%s" % (sampleDir, item['Content']['Path'])
            param['readtype'] = 'merged'
        elif name == 'Input.reads-sample':
            sample = item['Content']

            sampleDir = '/data/input/samples/%s/Data/Intensities/BaseCalls' % sample['Id']
            if not os.path.exists(sampleDir):
                sampleDir = '/data/input/samples/%s' % sample['Id']

            for root, dirs, files in os.walk(sampleDir):
                R1files = fnmatch.filter(files, '*_R1_*')
                R2files = fnmatch.filter(files, '*_R2_*')

            assert(len(R1files) == len(R2files))
            if len(R1files) != 1:
                print >> sys.stderr, "WARNING: Samples with multiple files are not supported yet. Only the first pair will be used"

            param['readtype'] = 'paired'
            param['reads'] = [os.path.join(sampleDir, R1files[0]),
                              os.path.join(sampleDir, R2files[0])]
        elif name == 'Input.tau':
            param['tau'] = safe_cast(item['Content'], int, 4)
        elif name == 'Input.loci':
            param['loci'] = item['Content']
        elif name == 'Input.organism':
            param['organism'] = item['Content']
        elif name == 'Input.app-session-name':
            param['AppSessionName'] = item['Content']
        elif name == 'Input.min-cluster-size':
            param["min-cluster-size"] = safe_cast(item['Content'], int, 5)
        elif name == 'Input.min-fillin':
            param["min-fillin"] = safe_cast(item['Content'], float, 0.6)
        elif name == 'Input.additionalfiles':
            param["additionalfiles"] = "save" in item['Items']
        elif name == 'Input.pseudogenes':
            param["pseudogenes"] = "pseudogenes" in item['Items']
        elif name == 'Input.Projects':
            projectID = item['Items'][0]['Id']

    param['outdir'] = '/data/output/appresults/%s/%s' % (projectID, "IgReC")
    mkdir_p(param['outdir'])

    igrec_path = os.path.join(igrec_dir, "igrec.py")

    if param['readtype'] == "merged":
        command = '%s -s "%s" ' % (igrec_path, param["reads"])
    elif param['readtype'] == "paired":
        command = '%s -1 "%s" -2 "%s" ' % (igrec_path, param["reads"][0], param["reads"][1])
    else:
        sys.exit(1)

    command += ' --tau=%(tau)d --loci=%(loci)s -o "%(outdir)s" --threads=%(threads)d --min-cluster-size=%(min-cluster-size)d --min-fillin=%(min-fillin)f' % param
    if not param["pseudogenes"]:
        command += " --no-pseudogenes"
    if param["additionalfiles"]:
        command += " --debug"

    # print "Command line: %s" % command
    os.system(command)

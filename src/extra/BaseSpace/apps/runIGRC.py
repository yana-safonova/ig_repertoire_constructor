#!/usr/bin/env python2

import os
import fnmatch
import sys
import json


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
    # load json file
    jsonfile = open('/data/input/AppSession.json')
    jsonObject = json.load(jsonfile)

    runID = jsonObject['Id']
    runName = jsonObject['Name']
    runName = runName.replace(' ', '_')

    # print jsonObject
    # showDir("/data")


    param = {}
    param["threads"] = 32
    param["min-size"] = 5
    param["tau"] = 4
    param["chain"] = "all"
    param["additionalfiles"] = True
    param["pseudogenes"] = True


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
                R1files = fnmatch.filter(files,'*_R1_*')
                R2files = fnmatch.filter(files,'*_R2_*')

            assert(len(R1files) == len(R2files))
            if len(R1files) != 1:
                print "WARNING: Samples with multiple files are not supported yet. Onlythe first file pair will be used"

            param['readtype'] = 'paired'
            param['reads'] = [os.path.join(sampleDir, R1files[0]),
                              os.path.join(sampleDir, R2files[0])]
        elif name == 'Input.tau':
            param['tau'] = int(item['Content'])
        elif name == 'Input.chain':
            param['chain'] = item['Content']
        elif name == 'Input.min-size':
            param["min-size"] = safe_cast(item['Content'], int, 5)
        elif name == 'Input.additionalfiles':
            param["additionalfiles"] = "save" in item['Items']
        elif name == 'Input.pseudogenes':
            param["pseudogenes"] = "pseudogenes" in item['Items']
        elif name == 'Input.Projects':
            projectID = item['Items'][0]['Id']

    param['outdir'] = '/data/output/appresults/%s/%s' % (projectID, runID)
    os.system('mkdir -p "%s"' % param['outdir'])
    print param

    igrc_dir = "/ig_repertoire_constructor"
    igrc_path = os.path.join(igrc_dir, "ig_repertoire_constructor.py")

    if param['readtype'] == "merged":
        command = "%s -s %s " % (igrc_path, param["reads"])
    elif param['readtype'] == "paired":
        command = "%s -1 %s -2 %s " % (igrc_path, param["reads"][0], param["reads"][1])
    else:
        sys.exit(1)

    command += " --tau=%(tau)d --chain=%(chain)s -o %(outdir)s --threads=%(threads)d --min-size=%(min-size)d" % param
    if not param["pseudogenes"]:
        command += " --no-pseudogenes"
    if param["additionalfiles"]:
        command += " --debug"

    print "Command line: %s" % command
    # os.system('cd %s; %s' % (igrc_dir, command))
    os.system(command)

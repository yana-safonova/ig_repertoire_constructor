#!/usr/bin/env python2

import os
import fnmatch
import sys
import json
import urllib


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

if __name__ == "__main__":
    #app specific definitions (not needed for personal app)
    parameter_list = []
    # load json file
    jsonfile = open('/data/input/AppSession.json')
    jsonObject = json.load(jsonfile)

    print jsonObject

    # determine the number of properties
    numberOfPropertyItems = len(jsonObject['Properties']['Items'])
    # loop over properties
    sampleID = []
    sampleHref = []
    sampleName = []
    sampleDir = []

    param = {}
    for item in jsonObject['Properties']['Items']:
        name = item['Name']

        if name == 'Input.input-reads-id':
            ParentAppResult = item['Content']['ParentAppResult']
            sampleDir = "/data/appresults/%s" % ParentAppResult['id'],

            param['reads'] = "%s/%s" % (sampleDir, item['Path'])
            param['readtype'] = 'merged'
        elif name == 'Input.input-reads-sample':
            assert(len(item['Items']) == 1)

            for sample in item['Items']:
                sampleID.append(sample['Id'])
                sampleHref.append(sample['Href'])
                sampleName.append(sample['Name'])
                sampleDir = '/data/input/samples/%s/Data/Intensities/BaseCalls' % sample['id']
                if not os.path.exists(sampleDir):
                    sampleDir = '/data/input/samples/%s' % sample['id']

                for root, dirs, files in os.walk(sampleDir):
                    R1files = fnmatch.filter(files,'*_R1_*')
                    R2files = fnmatch.filter(files,'*_R2_*')

                assert(len(R1files) == len(R2files) == 1)
            param['readtype'] = 'merged'
            param['reads'] = [R1files[0], R2files[0]]
        elif name == 'Input.input-tau-id':
            param['tau'] = int(item['Content'])
        elif name == 'Input.input-chain-id':
            param['chain'] = item['Content']
        elif name == 'Input.Projects':
            projectID = item['Items'][0]['Id']

    param['outdir'] = '/data/output/appresults/%s' % projectID
    os.system('mkdir -p "%s"' % param['outdir'])
    print param

    showDir("/data")


    igrc_path = "/ig_repertoire_constructor/ig_repertoire_constructor.py"

    if param['readtype'] == "merged":
        command = "%s -s %s " % (igrc_path, param["reads"])
    else:
        command = "%s -1 %s -2 %s " % (igrc_path, param["reads"][0], param["reads"][1])

    command += "--tau=%(out)d --chain=%(chain)s -o %(outdir)s" % param

    print "Command line: %s" % command
    os.system(command)

"""
    for index in range(numberOfPropertyItems):
    # set sample parameters
        if jsonObject['Properties']['Items'][index]['Name'] == 'Input.Samples':
            for sample in range(len(jsonObject['Properties']['Items'][index]['Items'])):
                sampleID.append(jsonObject['Properties']['Items'][index]['Items'][sample]['Id'])
                sampleHref.append(jsonObject['Properties']['Items'][index]['Items'][sample]['Href'])
                sampleName.append(jsonObject['Properties']['Items'][index]['Items'][sample]['Name'])
                sampleDir = '/data/input/samples/%s/Data/Intensities/BaseCalls' %(sampleID[sample])
                if not os.path.exists(sampleDir):
                    sampleDir = '/data/input/samples/%s' %(sampleID[sample])
                for root, dirs, files in os.walk(sampleDir[sample]):
                    R1files = fnmatch.filter(files,'*_R1_*')
                    R2files = fnmatch.filter(files,'*_R2_*')
                if len(R1files) != len(R2files):
                    print "number of R1 and R2 files do not match"
                    sys.exit()
                sampleOutDir = '/data/output/appresults/%s/%s' %(projectID,sampleName[sample])
                os.system('mkdir -p "%s"' %(sampleOutDir))
    # create output file and print parameters to output file (this is where you would run the command)
                file = '%s/parameters.csv' %(sampleOutDir)
                outFile = open(file ,'w')
                count = 0
                for parameter in parameter_list:
                    count += 1
                    outFile.write('%s,%s\n' %(count,parameter))
                outFile.close()
    #create metadata file for each appresult
                metadataObject = metadatajson()
                metaJsonObject = json.loads(metadataObject)
    #modify metadataObject
                metaJsonObject['Name'] = jsonObject['Properties']['Items'][index]['Items'][sample]['Id']
                metaJsonObject['Description'] = 'Sample Description'
                metaJsonObject['HrefAppSession'] = jsonObject['Href']
                for href in sampleHref:
                    metaJsonObject['Properties'][0]['Items'].append(href)

                metadataFile = '%s/_metadata.json' %(sampleOutDir)
                outMetadataFile = open(metadataFile, 'w')
                json.dump(metaJsonObject,outMetadataFile)
"""

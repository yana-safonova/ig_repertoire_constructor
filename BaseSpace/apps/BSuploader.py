#!/usr/bin/env python2

from BaseSpacePy.api.BaseSpaceAPI import BaseSpaceAPI
from argparse import ArgumentParser
import os
import os.path


if __name__ == "__main__":
    parser = ArgumentParser("BaseSpace data uploader")
    parser.add_argument("file",
                        type=str,
                        help="Source file for upload")
    parser.add_argument("--project", "-p",
                        type=str,
                        default="Upload",
                        help="Target project name (default: %(default)s)")
    parser.add_argument("--dir", "-d",
                        type=str,
                        default="upload",
                        help="Target upload directory (default: %(default)s)")
    parser.add_argument("--sdir", "-s",
                        type=str,
                        default="",
                        help="Target subdirectory (default: %(default)s)")
    parser.add_argument("--name", "-n",
                        type=str,
                        help="Target file name")
    args = parser.parse_args()


    print "Connecting to BaseSpace API server..."
    myAPI = BaseSpaceAPI()
    user = myAPI.getUserById('current')

    print "Connected to server, username: %s" % user

    p = myAPI.createProject(args.project)

    appResults = p.createAppResult(myAPI,
                                   args.dir,
                                   "",
                                   appSessionId="")

    myAppSession = appResults.AppSession
    myAppSession.setStatus(myAPI, "running", "uploading data")

    if not args.name or len(args.file):
        args.name = os.path.basename(args.file)
    print "Uploading file %(file)s to %(project)s:%(dir)s/%(sdir)s/%(name)s" % args.__dict__

    appResults.uploadFile(myAPI,
                          args.file,
                          args.name,
                          args.sdir,
                          "text/plain")

    myAppSession.setStatus(myAPI, "complete", "data uploaded")
    print "File successfully uploaded"

#!/bin/bash

REPO=eodus/mydock
REPO=docker.illumina.com/ig_repertoire_constructor/igrc

docker run -it -v /:/mnt ${REPO} bash -c "cd /ig_repertoire_constructor; git pull; ./prepare_cfg; make"


ID=`docker ps -l -q`
docker commit $ID ${REPO}

docker push ${REPO}



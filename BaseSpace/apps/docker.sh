#!/bin/bash

docker run -it -v /:/mnt eodus/mydock bash -c "cd /ig_repertoire_constructor; git pull"


ID=`docker ps -l -q`
docker commit $ID eodus/mydock

docker push eodus/mydock



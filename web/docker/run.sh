#!/bin/bash

set -euxo pipefail

docker build . -t igrec-web

RUN_DIR="/Bmo/ashlemov/igrec-web/runs"
UPLOAD_DIR="/Bmo/ashlemov/igrec-web/uploads"
REDIS_DIR="/Bmo/ashlemov/igrec-web/redis"
TMP_DIR="/Bmo/ashlemov/igrec-web/tmp"
DATA_DIR="/Nancy/data"
ROOT_DIR="/Bmo"
mkdir -p "${RUN_DIR}"
mkdir -p "${UPLOAD_DIR}"
mkdir -p "${TMP_DIR}"
mkdir -p "${REDIS_DIR}"

docker run -p 15284:8000 \
    --mount type=bind,source="${RUN_DIR}",destination=/opt/y-tools/web/runs \
    --mount type=bind,source="${TMP_DIR}",destination=/opt/y-tools/web/tmp \
    --mount type=bind,source="${UPLOAD_DIR}",destination=/opt/y-tools/web/uploads \
    --mount type=bind,source="${REDIS_DIR}",destination=/opt/y-tools/web/redis \
    --mount type=bind,source="${DATA_DIR}",destination=/data,readonly \
    --mount type=bind,source="${ROOT_DIR}",destination=/host-root,readonly \
    --mount type=bind,source="${UPLOAD_DIR}",destination=/uploads,readonly \
    igrec-web:latest

# --user $(id -u) \
# TODO Use --user option on running contatiner
# Use
# ssh -C -q -N -o ServerAliveInterval=60 <USERNAME>@chihua.cab.spbu.ru -L 8000:localhost:15284
# To connect to server


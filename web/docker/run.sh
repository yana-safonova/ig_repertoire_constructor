#!/bin/bash

set -euxo pipefail

docker build . -t igrec-web

RUN_DIR="/Bmo/ashlemov/igrec-web/runs"
UPLOAD_DIR="/Bmo/ashlemov/igrec-web/uploads"
mkdir -p "${RUN_DIR}"
mkdir -p "${UPLOAD_DIR}"

docker run -p 15284:8000 \
    --mount type=bind,source="${RUN_DIR}",destination=/opt/y-tools/web/runs \
    --mount type=bind,source="${UPLOAD_DIR}",destination=/opt/y-tools/web/tmp \
    --mount type=bind,source="${UPLOAD_DIR}",destination=/uploads \
    --mount type=bind,source=/,destination=/host-root,readonly \
    igrec-web:latest

# Use
# ssh -C -q -N -o ServerAliveInterval=60 <USERNAME>@chihua.cab.spbu.ru -L 8000:localhost:15284
# To connect to server


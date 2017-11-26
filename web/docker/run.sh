docker build . -t igrec-web

mkdir -p /Bmo/ashlemov/igrec-web/runs
mkdir -p /Bmo/ashlemov/igrec-web/uploads

docker run -p 15284:8000 --mount type=bind,source=/Bmo/ashlemov/igrec-web/runs,destination=/opt/y-tools/web/runs \
    --mount type=bind,source=/Bmo/ashlemov/igrec-web/uploads,destination=/opt/y-tools/web/tmp \
    --mount type=bind,source=/,destination=/host-root,readonly \
    igrec-web:latest

# Use
# ssh -C -q -N -o ServerAliveInterval=60 ashlemov@chihua.cab.spbu.ru -L 8000:localhost:15284
# To connect to server


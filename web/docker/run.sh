docker build . -t igrec-web

mkdir -p /Bmo/ashlemov/igrec-web/runs
mkdir -p /Bmo/ashlemov/igrec-web/uploads

docker run -p 15284:8000 --mount type=bind,source=/Bmo/ashlemov/igrec-web/runs,destination=/opt/y-tools/web/runs \
    --mount type=bind,source=/Bmo/ashlemov/igrec-web/uploads,destination=/opt/y-tools/web/tmp \
    --mount type=bind,source=/,destination=/host-root,readonly \
    igrec-web:latest


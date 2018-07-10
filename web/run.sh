#!/bin/bash

# https://stackoverflow.com/questions/18279540/bash-shutdown-hook-or-kill-all-background-processes-when-main-process-is-kille
trap 'jobs -p | xargs kill' EXIT

cd redis
redis-server &
cd ..

celery -A igrec-web.celery worker --loglevel info --statedb=worker.state --events &
# --events is necessary for task listing (aka inspect())
# To read about statedb: http://docs.celeryproject.org/en/latest/userguide/workers.html#persistent-revokes


# gunicorn --workers=8 igrec-web:app --worker-class eventlet &
gunicorn --workers=8 igrec-web:app --bind=:8000 --worker-class eventlet &

while pgrep -P "$BASHPID" > /dev/null; do
    wait
done

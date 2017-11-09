#!/bin/bash

# https://stackoverflow.com/questions/18279540/bash-shutdown-hook-or-kill-all-background-processes-when-main-process-is-kille
trap 'jobs -p | xargs kill' EXIT

redis-server &

celery -A igrec-web.celery worker &

gunicorn -w 8 igrec-web:app --worker-class eventlet &

while pgrep -P "$BASHPID" > /dev/null; do
    wait
done

#!/usr/bin/env bash
set -e

echo "building"
make -j8

echo "exporting artifacts"
python ./src/extra/serg_tools/copy_env.py env 1
echo $1 > env/GIT_REVISION
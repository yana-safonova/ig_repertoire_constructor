echo "building"
./prepare_cfg
make -j8

echo "exporting artifacts"
python ./src/extra/serg_tools/copy_env.py env 1
echo %build.vcs.number% > env/GIT_REVISION
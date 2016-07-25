#!/bin/bash
cd ../../;
# ./prepare_cfg
make vdj;
# make check_essential
build/release/bin/vdj_labeler_0.3 configs/vdj_labeler_0.3/configs.info;

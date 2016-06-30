#!/bin/bash
cd ../../; ./prepare_cfg
make vdj;
build/release/bin/vdj_labeler_0.3 configs/vdj_labeler/configs.info;

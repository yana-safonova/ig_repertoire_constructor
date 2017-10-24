#!/bin/bash
cd ../../; 
make shm;
./build/release/bin/shm_kmer_model configs/shm_kmer_model/config.info;

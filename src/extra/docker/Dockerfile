FROM ubuntu:latest

MAINTAINER Alexander Shlemov version: 0.1

RUN apt-get update && apt-get install -y git gcc g++ cmake python-matplotlib python-biopython && apt-get clean && rm -rf /var/lib/apt/lists/*
WORKDIR /opt
ENV upd 1
RUN git clone --depth=1 https://github.com/yana-safonova/ig_repertoire_constructor.git -b master
WORKDIR /opt/ig_repertoire_constructor


RUN ./prepare_cfg
RUN make

RUN ./igrec.py --test

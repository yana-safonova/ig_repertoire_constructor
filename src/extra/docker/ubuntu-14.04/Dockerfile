FROM ubuntu:14.04

MAINTAINER Alexander Shlemov version: 0.1

RUN apt-get update && apt-get install -y git gcc g++ cmake python-matplotlib python-biopython python-numpy python-scipy python-pandas python-pip
RUN python -m pip install --user seaborn
WORKDIR /opt
ENV upd 1
RUN git clone --depth=1 https://github.com/yana-safonova/ig_repertoire_constructor.git -b testing
WORKDIR /opt/ig_repertoire_constructor


RUN make -j3
RUN make check -j3

RUN make deb

RUN dpkg -i packages/*.deb

WORKDIR /opt

RUN igrec.py --test
RUN igquast.py --test --reference-free
RUN barcoded_igrec.py --test
RUN dense_subgraph_finder.py --test
RUN diversity_analyzer.py --test
RUN mass_spectra_analyzer.py --test

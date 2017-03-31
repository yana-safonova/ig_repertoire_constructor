FROM centos:7

MAINTAINER Alexander Shlemov version: 0.1

RUN yum -y install epel-release
RUN yum -y update
RUN yum -y install git gcc gcc-c++ make cmake python-biopython python-numpy python-pandas python-pip

RUN python -m pip install --user seaborn
# matplotlib and scipy from centos repo are obsolete
RUN python -m pip install --user --upgrade matplotlib
RUN python -m pip install --user --upgrade scipy

WORKDIR /opt
ENV upd 1
RUN git clone --depth=1 https://github.com/yana-safonova/ig_repertoire_constructor.git -b testing
WORKDIR /opt/ig_repertoire_constructor

RUN make -j3
RUN make check -j3

RUN yum -y install rpm-build
RUN make rpm

RUN rpm -i packages/*.rpm

WORKDIR /opt

RUN igrec.py --test
RUN igquast.py --test --reference-free
RUN barcoded_igrec.py --test
RUN dense_subgraph_finder.py --test
RUN diversity_analyzer.py --test
RUN mass_spectra_analyzer.py --test

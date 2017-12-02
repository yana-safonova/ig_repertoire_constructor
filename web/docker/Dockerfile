# Each instruction in this file generates a new layer that gets pushed to your local image cache
#

#
# Lines preceeded by # are regarded as comments and ignored
#

#
# The line below states we will base our new image on the Latest Official Ubuntu
FROM ubuntu:latest

#
# Identify the maintainer of an image
MAINTAINER Alexander Shlemov "shlemovalex@gmail.com"

#
# Update the image to the latest packages
RUN apt-get update && apt-get upgrade -y

#
# Install NGINX to test.
RUN apt-get install nginx -y

# Install python and pip
RUN apt-get install python-pip python-dev -y


RUN apt-get install wget -y
# Install Redis (taken from https://github.com/dockerfile/redis/blob/master/Dockerfile)
RUN \
  cd /tmp && \
  wget http://download.redis.io/redis-stable.tar.gz && \
  tar xvzf redis-stable.tar.gz && \
  cd redis-stable && \
  make && \
  make install && \
  cp -f src/redis-sentinel /usr/local/bin && \
  mkdir -p /etc/redis && \
  cp -f *.conf /etc/redis && \
  rm -rf /tmp/redis-stable* && \
  sed -i 's/^\(bind .*\)$/# \1/' /etc/redis/redis.conf && \
  sed -i 's/^\(daemonize .*\)$/# \1/' /etc/redis/redis.conf && \
  sed -i 's/^\(dir .*\)$/# \1\ndir \/data/' /etc/redis/redis.conf && \
  sed -i 's/^\(logfile .*\)$/# \1/' /etc/redis/redis.conf


RUN pip install --upgrade pip
RUN pip install flask flask_session flask_socketio flask_autoindex flask_silk celery redis gunicorn eventlet


# Install igrec

RUN apt-get install -y git gcc g++ cmake python-matplotlib python-biopython python-numpy python-scipy python-pandas python-pip
RUN python -m pip install --user seaborn
WORKDIR /opt
ENV upd 1
RUN git clone --depth=1 https://github.com/ablab/y-tools.git -b igrec_web
WORKDIR /opt/y-tools

RUN make -j3

# RUN igrec.py --test
# RUN igquast.py --test --reference-free
# RUN barcoded_igrec.py --test
# RUN dense_subgraph_finder.py --test
# RUN diversity_analyzer.py --test
# RUN mass_spectra_analyzer.py --test



#
# REMOVE
EXPOSE 80
EXPOSE 8080
EXPOSE 5000



# Run redis
# CMD ["redis-server", "/etc/redis/redis.conf"]
# Last is the actual command to start up NGINX within our Container
# CMD ["nginx", "-g", "daemon off;"]

WORKDIR /opt/y-tools/web
RUN git fetch origin webserver
RUN git fetch origin
RUN git fetch

RUN git config remote.origin.fetch "+refs/heads/*:refs/remotes/origin/*"
RUN git fetch origin webserver

RUN git checkout webserver

# We have to use special version of eventlet, since the latest one is buggy
# https://github.com/nameko/nameko/issues/425
RUN pip install eventlet==0.20.0

ENV VERSION 11
RUN git pull origin webserver

RUN mkdir -p /host-root

EXPOSE 8000
RUN rm /etc/redis/redis.conf
CMD ["bash", "run.sh"]





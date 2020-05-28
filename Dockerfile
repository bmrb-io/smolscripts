FROM frolvlad/alpine-glibc

MAINTAINER dmaziuk
EXPOSE 9090
WORKDIR /opt/wsgi

RUN apk update && \
    apk --no-cache add bash ca-certificates wget libxext libxrender libstdc++ && \
    update-ca-certificates && \
    apk --update add tzdata && \
    cp /usr/share/zoneinfo/America/Chicago /etc/localtime && \
    apk del tzdata

RUN  apk --no-cache add mc git perl-git git-svn

RUN echo 'export PATH=/opt/anaconda/bin:$PATH' > /etc/profile.d/anaconda.sh && \
    wget --quiet https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/anaconda.sh && \
    /bin/bash ~/anaconda.sh -b -p /opt/anaconda && \
    rm ~/anaconda.sh

ENV PATH /opt/anaconda/bin:$PATH

RUN conda install -y -c rdkit rdkit

RUN conda install -y -c openbabel openbabel

RUN conda install -y -c conda-forge libiconv uwsgi werkzeug

RUN conda clean --all --yes

RUN addgroup -g 101 -S uwsgi && \
    adduser -u 100 -S uwsgi -G uwsgi

RUN cd /opt/wsgi && \
    git svn clone http://octopus.bmrb.wisc.edu/svn/metabolomics_pipe . && \
    chown -R uwsgi:uwsgi .

CMD [ "uwsgi", \
    "--http", ":9090", \
    "--uid", "uwsgi", \
    "--gid", "uwsgi", \
    "--pidfile", "/run/uwsgi.pid", \
    "--python-autoreload", "10", \
    "--wsgi-file", "/opt/wsgi/wsgi.py" ]

#FROM ubuntu:20.04
FROM python:3.12-slim-bullseye

LABEL maintainer="Tobias Krause (krauset@mskcc.org), Hussen Mohammed Ibrahim (ibrahih3@mskcc.org)"

# install
COPY requirements.txt /requirements.txt
RUN pip install -r /requirements.txt

# copy
COPY src/* /opt/

# data
COPY ./data/3M-february-2018.txt.gz /opt/data/
COPY ./data/3M-3pgex-may-2023.txt.gz /opt/data/

# setup
WORKDIR /opt
CMD ["/bin/bash"]
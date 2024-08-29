#FROM ubuntu:20.04
FROM python:3.12-slim-bullseye

LABEL maintainer="Tobias Krause (krauset@mskcc.org)"

# install
COPY requirements.txt /requirements.txt
RUN pip install -r /requirements.txt

# copy
COPY src/* /opt/

# data
COPY ./data/3M-february-2018.txt.gz /opt/data/
COPY ./data/10x-hto-gex-mapper.pickle /opt/data/10x-hto-gex-mapper.pickle

# setup
WORKDIR /opt
CMD ["/bin/bash"]
FROM python:3.12-slim-bullseye

LABEL maintainer="Tobias Krause (krauset@mskcc.org)"

ARG VERSION_HTO_DND

# install stable version from pypi
RUN pip install hto==${VERSION_HTO_DND} --no-cache-dir

WORKDIR /opt

ENTRYPOINT [ "hto" ]

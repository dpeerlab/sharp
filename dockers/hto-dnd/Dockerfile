FROM python:3.12-slim-bullseye

LABEL maintainer="Tobias Krause (krauset@mskcc.org)"

WORKDIR /opt

# requirements
COPY requirements.txt /tmp/requirements.txt
RUN pip install --no-cache-dir -r /tmp/requirements.txt && \
    rm -rf /root/.cache/pip /tmp/requirements.txt

# hto-dnd
ARG VERSION_HTO_DND
RUN pip install --no-cache-dir hto==${VERSION_HTO_DND}

# Default entrypoint
ENTRYPOINT ["hto"]

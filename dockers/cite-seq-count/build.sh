#!/bin/bash -e

source config.sh

docker build \
    --tag ${image_name}:${version} \
    --platform linux/amd64 \
    --build-arg CITE_SEQ_COUNT_VERSION=${version} \
    -f Dockerfile-${version} .

#!/bin/bash

source config.sh

# build ${image_name}:${version}
docker build \
    --tag ${image_name}:${version} \
    --tag ${image_name}:latest \
    --platform linux/amd64 \
    --build-arg VERSION_HTO_DND=${version} .

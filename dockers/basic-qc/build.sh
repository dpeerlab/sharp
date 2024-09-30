#!/bin/bash -e

source config.sh
docker build --platform linux/amd64 -t ${image_name}:${version} .

#!/bin/bash

source config.sh

docker run --rm -it \
  -v $(pwd):/opt \
  -w /opt \
  ${image_name}:${version} \
  pytest -v /opt/tests \

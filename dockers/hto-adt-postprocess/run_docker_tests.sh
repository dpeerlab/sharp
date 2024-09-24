#!/bin/bash

source config.sh

docker run --rm -it \
  -v $(pwd):/opt/project \
  -w /opt/project \
  ${image_name}:${version} \
  pytest -v /opt/project/tests
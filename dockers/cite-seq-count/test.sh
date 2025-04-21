#!/bin/bash

source config.sh

docker run -it --rm ${image_name}:${version} CITE-seq-Count --help
docker run -it --rm ${image_name}:${version} CITE-seq-Count --version

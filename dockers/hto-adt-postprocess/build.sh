#!/bin/bash -e

source config.sh

if [ ! -r ./data/10x-hto-gex-mapper-3-prime.pickle ] || [ ! -r ./data/10x-hto-gex-mapper-5-prime.pickle ];
then
    # pre-build 10x-hto-gex-mapper.pickle
    python hto_gex_mapper.py
fi

docker build --platform linux/amd64 -t ${image_name}:${version} .

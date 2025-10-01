#!/bin/bash -e

source config.sh

echo "${registry}/${image_name}:${version}"

docker push ${registry}/${image_name}:${version}

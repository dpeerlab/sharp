#!/bin/bash -e

source config.sh

echo "${registry}/${image_cellranger_name}:${version}"

docker push ${registry}/${image_cellranger_name}:${version}
docker push ${registry}/${image_cromwell_name}:${version}

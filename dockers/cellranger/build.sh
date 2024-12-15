#!/bin/bash

source config.sh

# build ${image_name}:${version}
docker build \
    --tag ${image_cellranger_name}:${version} \
    --platform linux/amd64 \
    --build-arg DOWNLOAD_URL=${download_url} \
    --build-arg CELLRANGER_VERSION=${version} .

# hack: comment the ENTRYPOINT and CMD lines to make it work for cromwell
# this will generate `Dockerfile.cromwell` and build it under the name `cromwell-${image_name}:${version}`
# https://github.com/broadinstitute/cromwell/issues/2461
cat Dockerfile \
    | sed 's/^ENTRYPOINT \[/# ENTRYPOINT \[/g' \
    | sed 's/^CMD \[/# CMD \[/g' > Dockerfile.cromwell

# build cromwell-${image_name}:${version}
docker build \
    --tag ${image_cromwell_name}:${version} \
    --platform linux/amd64 \
    --build-arg DOWNLOAD_URL=${download_url} \
    --build-arg CELLRANGER_VERSION=${version} \
    -f Dockerfile.cromwell .

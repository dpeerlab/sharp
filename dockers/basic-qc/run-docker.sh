#!/bin/bash

source config.sh

sample_name="test_hashtag"
path_h5ad="data/hashtag/adata.h5ad"
path_report="data/hashtag/run_report.yaml"
path_reads="data/hashtag/reads"
path_outputs="data/outputs/hashtag"

docker run --rm \
    --platform linux/amd64 \
    --volume $(pwd)/data:/opt/data \
    ${image_name}:${version} \
    render hashtag_report.html \
            --sample-name $sample_name \
            --path-h5ad $path_h5ad \
            --path-report $path_report \
            --path-reads $path_reads \
            --path-output $path_outputs/hashtag_report.html

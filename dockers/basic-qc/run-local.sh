#!/bin/bash

source config.sh

# this cell is tagged `parameters` and will be overridden by papermill
sample_name="test_hashtag"
path_h5ad="data/hashtag/adata.h5ad"
path_report="data/hashtag/run_report.yaml"
path_reads="data/hashtag/reads"
path_outputs="outputs/hashtag"

mkdir -p ${path_outputs}

papermill \
    $path_notebook_hashtag ${path_outputs}/hashtag_report.ipynb \
    --parameters sample_name $sample_name \
    --parameters path_h5ad $path_h5ad \
    --parameters path_report $path_report \
    --parameters path_reads $path_reads \
    --stdout-file ${path_outputs}/hashtag.stdout.txt \
    --log-output

jupyter nbconvert ${path_outputs/hashtag_report.ipynb} \
    --to html \
    --sanitize-html \
    --theme=light \
    --output ${path_outputs}/hashtag_report.html

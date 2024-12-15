#!/bin/bash

source config.sh

templateNotebook=notebooks/inspect-hashtag-v4.ipynb
sampleName=ZZ-000_cite_test
h5ad=test/hashtag-v4/adata_final.h5ad
path_reads=test/hashtag-v4/reads
runReport=test/hashtag-v4/run_report.yaml
path_outdir=test/hashtag-v4/outputs

docker run --rm \
    -p 8888:8888 \
    -v $(pwd):/inputs \
    ${image_name}:${version} \
    /bin/bash -c " \
    cd /inputs && \
    papermill \
        $templateNotebook $path_outdir/$sampleName.QC.ipynb \
        --parameters sample_name $sampleName \
        --parameters path_h5ad $h5ad \
        --parameters path_report $runReport \
        --parameters path_reads $path_reads \
        --parameters path_outdir $path_outdir \
        --stdout-file $path_outdir/$sampleName.QC.stdout.txt \
        --log-output && \
    jupyter nbconvert --to html_toc --ExtractOutputPreprocessor.enabled=False ${path_outdir}/${sampleName}.QC.ipynb \
    "

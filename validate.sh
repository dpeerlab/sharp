#!/usr/bin/env bash

echo ">> Hashtag 10x v3 / TotalSeq-A"
java -jar ~/Applications/womtool.jar \
    validate \
    Hashtag.wdl \
    --inputs ./configs/hashtag-10x-v3-tsa.inputs.json

echo ">> Hashtag 10x v3 / TotalSeq-B barcode translation"
java -jar ~/Applications/womtool.jar \
    validate \
    Hashtag.wdl \
    --inputs ./configs/hashtag-10x-v3-tsb.inputs.json

echo ">> Hashtag InDrop Methanol"
java -jar ~/Applications/womtool.jar \
    validate \
    Hashtag.wdl \
    --inputs ./configs/hashtag-indrop-methanol.inputs.json

echo ">> CITE-seq"
java -jar ~/Applications/womtool.jar \
    validate \
    CiteSeq.wdl \
    --inputs ./configs/citeseq.inputs.json

echo ">> CellPlex"
java -jar ~/Applications/womtool.jar \
    validate \
    Hashtag.wdl \
    --inputs ./configs/cellplex.inputs.json

echo ">> ASAP-seq"
java -jar ~/Applications/womtool.jar \
    validate \
    AsapSeq.wdl \
    --inputs ./configs/asapseq-tsa.inputs.json

echo ">> Hashtag 10x v2 / 5'"
java -jar ~/Applications/womtool.jar \
    validate \
    Hashtag.wdl \
    --inputs ./configs/hashtag-10x-v2-5p.inputs.json

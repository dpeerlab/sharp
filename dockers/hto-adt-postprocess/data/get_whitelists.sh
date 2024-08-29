cd /Users/krauset/projects/sharp/dockers/hto-adt-postprocess/0.3.8

# setup data
# laucnh docker
docker run -d -t sailmskcc/cromwell-cellranger:8.0.0

# prep data
rm -r data/*.txt
rm -r data/*.txt.gz

# copy data
docker cp $(docker ps -lq):/opt/cellranger-8.0.0/lib/python/cellranger/barcodes/translation/3M-3pgex-may-2023.txt.gz data
docker cp $(docker ps -lq):/opt/cellranger-8.0.0/lib/python/cellranger/barcodes/translation/3M-february-2018.txt.gz data

# copy barcodes
cp barcodes/3M-3pgex-may-2023.txt.gz ../data
cp barcodes/3M-5pgex-jan-2023.txt.gz ../data 
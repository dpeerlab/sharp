PATH_CR=/Users/krauset/projects/sharp_dev/dockers/cellranger

# get download-url from
# https://www.10xgenomics.com/support/software/cell-ranger/downloads

# download to /tmp
cd /tmp
curl -o cellranger-9.0.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-9.0.0.tar.gz?Expires=1734263545&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=CFs--LhYp4b9U4ZCArF1ku1HkL-IEuBxICFu8OBo2GhxHPomNdKZS~GC-SmDy7s6Ahh7sND-wF1GE2gkll5JLI5jv7yvEw2fn99BwsI4X~zXyDWM~QfXv0uYVSLb2I8K~zrGmOi9tJpeqe~CZ4YHRgh7ojHBvMaDQK4Ldz8q0HEql79~Ph8A16pMBA9HbydbQHsNtBshte1A59RVqxBVtDE5rNunsOeBbYIPfN9Sbjw7335ASDA28KkAcvCfsy5ECilCGaxuubwVLyA0OTsDbJnMx~1CyhfxoLyVuneQMBGQXw~JZI62mEjZbmELJdAuTUcSu5rc2Lb~8EmwQJeVcw__"

# extract to $PATH_CR
tar -xzvf cellranger-9.0.0.tar.gz -C $PATH_CR
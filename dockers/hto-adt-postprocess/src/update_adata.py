#!/usr/bin/env python
# coding: utf-8

import sys
import argparse
import logging
from anndata._core.anndata import AnnData
import anndata as ad
import pandas as pd

numba_logger = logging.getLogger("numba")
numba_logger.setLevel(logging.WARNING)

logger = logging.getLogger("updata_adata")

logging.basicConfig(
    level=logging.DEBUG,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler("update_adata.log"),
        logging.StreamHandler(sys.stdout),
    ],
)

def updata_adata(
    path_class: str,
    path_adata_in: str,
    path_adata_out: str,
):

    logger.info(f"Loading AnnData {path_adata_in}...")
    adata = ad.read_h5ad(path_adata_in)

    logger.info("Loading classification...")
    df_class = pd.read_csv(path_class, sep="\t", index_col=0, compression="gzip")

    logger.info("Adding classification to AnnData...")
    adata.obs["hash_id"] = pd.Categorical(df_class.hashID)

    logger.info(f"Writing AnnData to {path_adata_out}...")
    adata.write(path_adata_out)


def parse_arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--class",
        action="store",
        dest="path_class",
        help="path to hashtag classification file (*.tsv.gz)",
        required=True,
    )

    parser.add_argument(
        "--adata-out",
        action="store",
        dest="path_adata_out",
        help="path to output AnnData (.h5ad)",
        required=True,
    )

    parser.add_argument(
        "--adata-in",
        action="store",
        dest="path_adata_in",
        help="path to input AnnData (.h5ad)",
        required=True,
    )

    # parse arguments
    params = parser.parse_args()

    return params


if __name__ == "__main__":

    params = parse_arguments()

    logger.info("Starting...")

    updata_adata(
        path_class=params.path_class,
        path_adata_in=params.path_adata_in,
        path_adata_out=params.path_adata_out,
    )

    logger.info("DONE.")

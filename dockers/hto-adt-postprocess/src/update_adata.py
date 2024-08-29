#!/usr/bin/env python
# coding: utf-8

import sys
import argparse
import logging
from anndata._core.anndata import AnnData
import anndata as ad
import pandas as pd
from dna3bit import DNA3Bit
from translate_barcodes import translate_barcodes

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


def translate(adata: AnnData, chemistry: str):

    # translate (TotalSeq-B/C HTO <--> GEX)
    translated_barcodes = translate_barcodes(
        adata.obs["barcode_sequence"].values, chemistry=chemistry
    )

    # encode nucleotide barcodes into numerical barcodes
    dna3bit = DNA3Bit()
    numerical_barcodes = list(
        map(lambda x: str(dna3bit.encode(x)), translated_barcodes)
    )
    adata.obs_names = numerical_barcodes


def updata_adata(
    path_class: str,
    path_adata_in: str,
    path_adata_out: str,
    translate_10x_barcodes: bool,
    chemistry: str,
):

    logger.info(f"Loading AnnData {path_adata_in}...")
    adata = ad.read_h5ad(path_adata_in)

    logger.info("Loading classification...")
    df_class = pd.read_csv(path_class, sep="\t", index_col=0, compression="gzip")

    logger.info("Adding classification to AnnData...")
    adata.obs["hash_id"] = pd.Categorical(df_class.hashID)

    if translate_10x_barcodes:
        logger.info("Translating TotalSeq-B/C HTO <--> GEX barcodes...")
        translate(adata, chemistry=chemistry)

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

    parser.add_argument(
        "--10x-barcode-translation",
        action="store_true",
        dest="translate_10x_barcodes",
        help="Translate HTO barcodes to GEX barcodes",
        default=False,
    )

    parser.add_argument(
        "--chemistry",
        action="store",
        dest="chemistry",
        help="Chemistry, as specified in the emulsion sheet, helps determine the whitelist.",
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
        translate_10x_barcodes=params.translate_10x_barcodes,
        chemistry=params.chemistry,
    )

    logger.info("DONE.")

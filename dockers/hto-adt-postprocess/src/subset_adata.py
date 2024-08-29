#!/usr/bin/env python
# coding: utf-8

import sys
import argparse
import logging
import anndata as ad
import pandas as pd

from dna3bit import DNA3Bit

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

def assert_cellbarcodes(barcodes, top_k=1000):
    """
    Assert that cell barcodes only contain ACGT and unique.
    """
    for i, barcode in enumerate(barcodes[:top_k]):
        if not barcode.isalnum():
            logger.warning(f"Cell barcode '{barcode}' contains non-alphanumeric characters.")
    assert len(barcodes) == len(set(barcodes)), "Cell barcodes are not unique."

def convert_barcodes(x):
    """
    Convert a cell barcode to DNA3Bit if all letters are numeric.
    """
    encoder_decoder = DNA3Bit()
    x = map(lambda i: encoder_decoder.decode(int(i)).decode(), x)
    return list(x)

def subset_adata(
    path_adata_in: str,
    path_adata_out: str,
    path_cb_whitelist: str,
    convert: bool = True,
):
    """
    Subset an AnnData object to only include cells in a given whitelist. This only works with alphanumeric cell barcodes.
    """
        
    
    logger.info(f"Loading AnnData {path_adata_in}...")
    adata = ad.read_h5ad(path_adata_in)

    logger.info(f"Loading cell barcode whitelist {path_cb_whitelist}...")
    cb_whitelist = pd.read_csv(path_cb_whitelist, header=None, index_col=0).index.values

    logger.info(f"Converting cell barcode whitelist to DNA3Bit...")
    if convert:
        adata.obs_names = convert_barcodes(adata.obs_names)

    logger.info("Asserting cell barcodes...")
    assert_cellbarcodes(adata.obs_names)
    assert_cellbarcodes(cb_whitelist)

    logger.info("Subsetting AnnData...")
    cb_whitelist_clean = list(set(cb_whitelist).intersection(set(adata.obs_names)))
    adata = adata[cb_whitelist_clean]
    
    difference = set(cb_whitelist) - set(cb_whitelist_clean)
    if len(difference) > 0:
        logger.warning(f"Cell barcodes ({len(difference)}) are not in the AnnData object:  {' '.join(difference)}. Ignoring them.")
        
    logger.info(f"Writing AnnData to {path_adata_out}...")
    adata.write(path_adata_out)

def parse_arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--adata-in",
        action="store",
        dest="path_adata_in",
        help="path to input AnnData (.h5ad)",
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
        "--cb-whitelist",
        action="store",
        dest="path_cb_whitelist",
        help="path to cell barcode whitelist",
        required=True,
    )

    parser.add_argument(
        "--convert",
        action="store_true",
        dest="convert",
        help="convert cell barcodes to DNA3Bit",
        required=False,
    )

    # parse arguments
    params = parser.parse_args()

    return params


if __name__ == "__main__":

    params = parse_arguments()

    logger.info("Starting...")

    subset_adata(
        path_adata_in=params.path_adata_in,
        path_adata_out=params.path_adata_out,
        path_cb_whitelist=params.path_cb_whitelist,
        convert=params.convert,
    )

    logger.info("DONE.")

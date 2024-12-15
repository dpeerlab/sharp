#!/usr/bin/env python
# coding: utf-8

import sys
import shutil
import argparse
import logging
import pandas as pd

logger = logging.getLogger("fix_barcodes")

logging.basicConfig(
    level=logging.DEBUG,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler("fix_barcodes.log"),
        logging.StreamHandler(sys.stdout),
    ],
)


def fix_barcodes(
    path_counts_in: str,
    path_counts_out: str,
    cb_whitelist_method: str,
):
    """
    Subset an AnnData object to only include cells in a given whitelist. This only works with alphanumeric cell barcodes.
    """

    logger.info("Copy data to output path...")
    shutil.copytree(path_counts_in, path_counts_out)

    logger.info("Preparing barcodes...")
    path_barcodes = f"{path_counts_out}/barcodes.tsv.gz"

    if cb_whitelist_method.lower() == "10x":
        logger.info(f"Adding '-1' to barcodes...")
        barcodes = pd.read_csv(path_barcodes, header=None, names=["barcode"])
        barcodes["barcode"] = barcodes["barcode"] + "-1"
        barcodes.to_csv(path_barcodes, header=False, index=False)


def parse_arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--counts-in",
        action="store",
        dest="path_counts_in",
        help="path to input counts matrix, containing barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz",
        required=True,
    )

    parser.add_argument(
        "--counts-out",
        action="store",
        dest="path_counts_out",
        help="path to output counts matrix, containing barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz",
        required=True,
    )

    parser.add_argument(
        "--cb-whitelist-method",
        action="store",
        dest="cb_whitelist_method",
        help="method of whitelist generation. If 10x, add -1 to barcodes.",
        required=True,
    )

    # parse arguments
    params = parser.parse_args()

    return params


if __name__ == "__main__":

    params = parse_arguments()

    logger.info("Starting...")

    fix_barcodes(
        path_counts_in=params.path_counts_in,
        path_counts_out=params.path_counts_out,
        cb_whitelist_method=params.cb_whitelist_method,
    )

    logger.info("DONE.")

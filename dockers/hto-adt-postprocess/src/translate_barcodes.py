#!/usr/bin/env python
import os
import sys
import argparse
import pandas as pd
import numpy as np
import json
import csv
import gzip
import logging
from tqdm import tqdm
import hto_gex_mapper

logger = logging.getLogger("translate_barcodes")

logging.basicConfig(
    level=logging.DEBUG,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler("translate_barcodes.log"),
        logging.StreamHandler(sys.stdout),
    ],
)


def convert(df, path_hto_gex_mapper):

    # load pre-built 10x HTO <--> GEX mapper
    mapper = hto_gex_mapper.load(path_hto_gex_mapper)

    # translate
    translated_barcodes = df.index.map(lambda x: mapper[x])

    df2 = df.copy()
    df2.index = translated_barcodes

    return df2


def translate(
    path_barcodes,
    chemistry,
):
    barcodes = pd.read_csv(
        path_barcodes, sep="\t", index_col=0, header=None, compression="gzip"
    )

    logger.info("Loaded barcodes ({})".format(len(barcodes)))

    # evaluate which whitelist to use
    logger.info(f"Determining which whitelist to use for {chemistry}...")
    path_hto_gex_mapper = hto_gex_mapper.decide_which_whitelist(chemistry)

    # translate HTO barcodes to GEX barcodes
    logger.info("Translating TotalSeq-B/C HTO <--> GEX barcodes...")
    df_final = convert(barcodes, path_hto_gex_mapper=path_hto_gex_mapper)

    df_final.to_csv("barcodes-translated.tsv.gz", header=None, compression="gzip")


def parse_arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--barcodes",
        action="store",
        dest="path_barcodes",
        help="path to barcode file (e.g. 10x's barcodes.tsv.gz)",
        required=True,
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

    translate(
        path_barcodes=params.path_barcodes,
        chemistry=params.chemistry,
    )

    logger.info("DONE.")

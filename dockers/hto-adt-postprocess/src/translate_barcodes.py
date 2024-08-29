#!/usr/bin/env python
import sys
import argparse
import pandas as pd
import logging
from hto_gex_mapper import decide_which_whitelist

logger = logging.getLogger("translate_barcodes")

logging.basicConfig(
    level=logging.DEBUG,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler("translate_barcodes.log"),
        logging.StreamHandler(sys.stdout),
    ],
)


def translate_barcodes(barcodes, chemistry: str):

    # get whitelist
    path_translation = decide_which_whitelist(chemistry)

    # translate
    translation_df = pd.read_csv(path_translation, sep="\t", index_col=0, header=None)

    translated_barcodes = translation_df.loc[barcodes].values.flatten()

    return translated_barcodes


def convert(df, chemistry: str):

    # translate
    index_new = translate_barcodes(df.index, chemistry)
    df_out = df.set_index(index_new)

    return df_out


def translate(
    path_barcodes,
    chemistry,
):
    barcodes = pd.read_csv(
        path_barcodes, sep="\t", index_col=0, header=None, compression="gzip"
    )

    logger.info("Loaded barcodes ({})".format(len(barcodes)))

    # translate HTO barcodes to GEX barcodes
    logger.info("Translating TotalSeq-B/C HTO <--> GEX barcodes...")
    df_final = convert(barcodes, chemistry)

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

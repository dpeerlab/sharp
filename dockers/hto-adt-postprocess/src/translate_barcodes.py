"""
Translate barcodes from HTO <--> GEX (whitelists are symmetrical)
"""

#!/usr/bin/env python
import sys
import argparse
import pandas as pd
import anndata as ad
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
    data_type="pandas",
    output_path=None
):
    if data_type == "pandas":
        barcodes = pd.read_csv(path_barcodes, sep="\t", index_col=0, header=None, compression="gzip")
        df_final = convert(barcodes, chemistry)
        if output_path is None:
            output_path = "barcodes-translated.tsv.gz"
        df_final.to_csv(output_path, header=None, compression="gzip")

    elif data_type == "adata":
        adata = ad.read_h5ad(path_barcodes)
        adata.obs = convert(adata.obs, chemistry)
        if output_path is None:
            output_path = "adata_translated.h5ad"
        adata.write_h5ad(output_path)

    else:
        raise ValueError(f"Currently supports only 'pandas' and 'adata'. Got '{data_type}'")


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

    parser.add_argument(
        "--data-type",
        action="store",
        dest="data_type",
        help="Type of data to translate. Currently supports 'pandas' and 'adata'.",
        required=False,
        default="pandas",
    )

    parser.add_argument(
        "--output-path",
        action="store",
        dest="output_path",
        help="Path to output file. If not specified, will be saved in the current directory.",
        required=False,
        default=None,
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
        data_type=params.data_type,
        output_path=params.output_path,
    )

    logger.info("DONE.")

"""
Translate barcodes from HTO <--> GEX (whitelists are symmetrical)
"""

#!/usr/bin/env python
import os
import sys
import argparse
import json
from typing import Union
import pandas as pd
import anndata as ad
import logging

logger = logging.getLogger("translate_barcodes")

logging.basicConfig(
    level=logging.DEBUG,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler("translate_barcodes.log"),
        logging.StreamHandler(sys.stdout),
    ],
)

def decide_which_whitelist(
    chemistry: str=None,
    whitelist_name: str=None,
    path_whitelists="/opt/data"
):
    """
    Based on 10x information, decide which whitelist to use. (https://kb.10xgenomics.com/hc/en-us/articles/115004506263-What-is-a-barcode-whitelist-)

    Args:
    - chemistry: str, chemistry used in the experiment. E.g. "V2", "V3", "V3.1", "V3.1.1". Uses the whitelist.json file to decide which whitelist to use.
    - whitelist_name: str, name of the whitelist file. If not specified, will be decided based on the chemistry.
    - path_whitelists: str, path to the whitelist json file. /opt/data is the default used in the docker container.
    """

    assert not (chemistry is None and whitelist_name is None), "Either chemistry or whitelist_name must be specified."
    assert chemistry is None or whitelist_name is None, "Only one of chemistry or whitelist_name must be specified."

    # get whitelist
    if chemistry is not None:
        # read
        with open(os.path.join(path_whitelists, "whitelists.json")) as fin:
            whitelists = json.load(fin)
        # get
        found = False
        for whitelist_name, chemistries in whitelists.items():
            if chemistry in chemistries:
                found = True
                break
        if not found:
            raise ValueError(f"Chemistry {chemistry} not found in the whitelist.json file.")
        path_whitelist = os.path.join(path_whitelists, whitelist_name)
    else:
        path_whitelist = os.path.join(path_whitelists, whitelist)

    assert os.path.exists(path_whitelist), f"Whitelist file '{whitelist_name}' does not exist in '{path_whitelists}'."
    return path_whitelist


def translate_barcodes(
    barcodes: Union[list, pd.Series, pd.Index],
    chemistry: str,
    base_path: str="/opt"
):
    """
    Transl
    """
    # to series
    allowed = [list, pd.Series, pd.Index]
    assert any([isinstance(barcodes, a) for a in allowed]), f"barcodes must be one of {allowed}. Got {type(barcodes)}"
    if isinstance(barcodes, list):
        barcodes = pd.Series(barcodes)

    # get whitelist
    path_translation = decide_which_whitelist(chemistry, path_whitelists=f"{base_path}/data")

    # remove -1 suffix
    barcodes = barcodes.str.replace("-1", "")

    # translate
    translation_df = pd.read_csv(path_translation, sep="\t", index_col=0, header=None)

    missing = list(set(barcodes) - set(translation_df.index))
    nmissing = len(missing)
    if nmissing > 0:
        raise ValueError(f"'{nmissing}' out of '{len(barcodes)}' barcodes are not in the whitelist. Got for example: {' '.join(missing[:5])}")
    translated_barcodes = translation_df.loc[barcodes].values.flatten()

    return translated_barcodes

def translate(
    path_barcodes,
    chemistry,
    data_type="pandas",
    output_path=None,
    base_path="/opt",
    kwargs_read_csv={},
):
    assert isinstance(path_barcodes, str), f"path_barcodes must be a string. Got {type(path_barcodes)}"
    if data_type == "pandas":
        # read args
        kwargs_read_csv["header"] = kwargs_read_csv.get("header", None)
        kwargs_read_csv["index_col"] = kwargs_read_csv.get("index_col", 0)
        kwargs_read_csv["sep"] = kwargs_read_csv.get("sep", "\t")
        # read
        barcodes = pd.read_csv(path_barcodes, **kwargs_read_csv)
        # translate
        barcodes.index = translate_barcodes(barcodes.index, chemistry, base_path=base_path)
        if output_path is None:
            output_path = "barcodes-translated.tsv.gz"
        barcodes.to_csv(output_path, header=None)

    elif data_type == "adata":
        # read or get adata
        adata = ad.read_h5ad(path_barcodes)
        adata.obs_names = translate_barcodes(adata.obs_names, chemistry, base_path=base_path)
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

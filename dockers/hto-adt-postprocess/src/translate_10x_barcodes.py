import os
import sys
import json
import argparse
import pandas as pd
import humanfriendly
import logging
from dna3bit import DNA3Bit


logger = logging.getLogger()

logging.basicConfig(
    level=logging.DEBUG,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler("translate_10x_barcodes.log"),
        logging.StreamHandler(sys.stdout),
    ],
)


def decode(barcodes):
    encoder_decoder = DNA3Bit()

    decoded = set(map(lambda x: encoder_decoder.decode(x).decode(), barcodes))

    return decoded


def decide_which_whitelist(chemistry, base_path="/opt"):
    """
    Based on 10x information, decide which whitelist to use. (https://kb.10xgenomics.com/hc/en-us/articles/115004506263-What-is-a-barcode-whitelist-)
    Chemistries are defined in the emulsion sheet.
    Not all whitelists are supported yet.

    Argument `base_path` describes where to find the whitelist json. Will be
    changed when tested locally with `pytest --local`.
    """
    with open(os.path.join(base_path, "data/whitelists.json")) as fin:
        whitelists = json.load(fin)

    for file, chemistries in whitelists.items():
        if chemistry in chemistries:
            return os.path.join(base_path, "data", file)
    else:
        raise ValueError("Chemistry {} not supported yet.".format(chemistry))


def translate(
    path_input, chemistry, separator, has_header, debug=False, base_path="/opt"
):

    df = pd.read_csv(path_input, sep=separator, header=0 if has_header else None)

    barcodes = df.iloc[:, 0].values

    if type(barcodes[0]) is str:
        # remove the suffix -1 in case Cell Ranger output
        barcodes = set(map(lambda x: x.strip("-1").strip(), barcodes))
    else:
        # SEQC outputs numerical barcode.
        # decode back to nucleotide sequence
        barcodes = decode(barcodes)

    # write translated barcodes
    path_10x_whitelist = decide_which_whitelist(chemistry, base_path=base_path)
    whitelist = pd.read_csv(path_10x_whitelist, sep="\t", header=None)
    whitelist = whitelist.loc[whitelist.iloc[:, 0].isin(barcodes)]  # subset
    whitelist.iloc[:, [1]].to_csv("translated-barcodes.txt", index=False, header=None)

    # outputs
    n = len(barcodes)
    m = len(whitelist)
    logger.info("Pre-translated: " + humanfriendly.format_number(n))
    logger.info("Post-translated: " + humanfriendly.format_number(m))

    if n != m:
        logger.error("Number of pre-/post-translated barcodes is different!")
        if debug:
            print(barcodes)
            print(whitelist)
            raise ValueError(
                f"Number of pre-/post-translated barcodes is different! {n} vs {m}"
            )
        else:
            exit(1)


def parse_arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--input-file",
        action="store",
        dest="path_input",
        help="path to barcode file (one per line) or multidimensional matrix (csv or tsv or gzipped)",
        required=True,
    )

    parser.add_argument(
        "--separator",
        action="store",
        dest="separator",
        help="a separator to be used when reading a matrix file",
        default=",",
    )

    parser.add_argument(
        "--header",
        action="store_true",
        dest="has_header",
        help="whether header exists",
        default=None,
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
        path_input=params.path_input,
        chemistry=params.chemistry,
        separator=params.separator,
        has_header=params.has_header,
    )

    logger.info("DONE.")

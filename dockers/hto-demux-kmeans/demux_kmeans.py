#!/usr/bin/env python
# coding: utf-8

import sys
import os
import argparse
import pandas as pd
import numpy as np
import yaml
import logging
import scipy.io
import scipy.stats
from sklearn.cluster import KMeans
from dna3bit import DNA3Bit
import warnings
from scipy.stats.mstats import gmean


logger = logging.getLogger("demux_kmeans")

logging.basicConfig(
    level=logging.DEBUG,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler("demux_kmeans.log"),
        logging.StreamHandler(sys.stdout),
    ],
)


def hto_demux(
    path_hto_umi_count_dir: str,
    mode: int,
    min_count_threshold: int
):
    # Construct sparse matrix of HTO counts
    matrix = scipy.io.mmread(
        os.path.join(path_hto_umi_count_dir, "matrix.mtx.gz")
        )
    barcodes = pd.read_csv(
        os.path.join(path_hto_umi_count_dir, "barcodes.tsv.gz"),
        header=None
    )[0]
    features = pd.read_csv(
        os.path.join(path_hto_umi_count_dir, "features.tsv.gz"),
        header=None
    )[0]

    # Remove barcodes with less than minimum counts
    negative_mask = np.ravel(matrix.sum(axis=0) > min_count_threshold)
    csr = matrix.tocsr()[:, negative_mask]

    # convert to numeric cell barcode
    dna3bit = DNA3Bit()
    numeric_barcodes = barcodes.apply(dna3bit.encode)

    # Convert to DataFrame
    df_umi = pd.DataFrame(
        csr.todense(),
        columns=numeric_barcodes[negative_mask],
        index=features,
    ).T

    logger.info(
        "Loaded HTO UMI count matrix %s",
        df_umi.shape[0]
    )

    # drop the column `unmapped`
    df_umi = df_umi.iloc[:, 0:-1]

    logger.info(f"Running in mode {mode}...")
    if mode == 1:
        # centered log-ratio (CLR) transformation
        df_clr = df_umi.apply(
            lambda row: np.log1p((row + 1) / gmean(row + 1)),
            axis=1
        )
    elif mode == 2:
        # very noisy methanol-based
        df_clr = df_umi.apply(lambda row: row - np.mean(row), axis=1)
        df_clr = df_clr.applymap(lambda x: 0 if x < 0 else x)
        df_clr = df_clr.apply(
            lambda row: np.log1p((row + 1) / gmean(row + 1)),
            axis=1
        )
    elif mode == 3:
        # aggresively rescue from doublets if in doubt
        df_clr = df_umi.apply(
            lambda row: row / gmean(row + 1),
            axis=1
        )
    else:
        raise Exception("Unrecognized mode...")

    # change column name to column index so that we can access by e.g. x[1]
    df_tmp = df_umi
    df_tmp.columns = range(0, len(df_tmp.columns))

    # for each row/barcode, get the index of the one with the largest UMI count
    ss_umi_largest = df_tmp.idxmax(axis=1)

    def kmeans_per_row(row):
        x = np.array(row).reshape(-1, 1)
        kmeans = KMeans(n_clusters=2, random_state=0).fit(x)
        y_predict = kmeans.predict(x)
        return y_predict

    logger.info("Running K-means...")
    # 227922838763364    [0, 1, 1, 1]
    # 239596337850148    [0, 0, 1, 0]
    # 164759051090203    [0, 1, 1, 1]
    # 191020391422693    [0, 1, 1, 0]
    # 204968413023541    [0, 0, 0, 1]
    with warnings.catch_warnings():
        # avoid ConvergenceWarning: Number of distinct clusters (1)
        # found smaller than n_clusters (2). 
        # Possibly due to duplicate points in X.
        warnings.simplefilter("ignore")
        df_kmeans = df_clr.apply(lambda row: kmeans_per_row(row), axis=1)

    df_kmeans_hotencoded = df_kmeans.apply(
        lambda x: "".join(str(y) for y in x)
    ).to_frame()

    # shorten and replace _ with -
    # ['HTO-301', 'HTO-302', 'HTO-303', 'HTO-304']
    hto_names = list(
        map(lambda name: name.split("-")[0].replace("_", "-"), df_clr.columns)
    )

    def demux_pass2(cb):

        # index of hto having the largest UMIs: 0, 1, 2, or 3
        idmax = ss_umi_largest[cb]

        # which group belongs to? 0 or 1
        group_id = df_kmeans_hotencoded.loc[cb][0][idmax]

        # how many hto belong that group?
        num_htos = df_kmeans_hotencoded.loc[cb][0].count(group_id)

        # if greater than or equal to two HTOs belong to that group,
        # it means doublet
        # return "Doublet" if doublet, return HTO ID if singlet
        # return "Doublet" if num_htos >= 2 else "Singlet"
        return "Doublet" if num_htos >= 2 else hto_names[idmax]

    df_class = pd.DataFrame(
        list(map(lambda cb: (cb, demux_pass2(cb)), df_kmeans_hotencoded.index))
    )
    df_class.columns = ["CB", "hashID"]
    df_class.set_index("CB", inplace=True)

    logger.debug(df_class.groupby(by="hashID").size())

    df_neg = pd.DataFrame(
        "Negative",
        index=numeric_barcodes[~negative_mask],
        columns=["hashID"],
    )
    df_neg.index.rename("CB", inplace=True)
    df_class = pd.concat([df_class, df_neg]).loc[numeric_barcodes]

    df_class.to_csv("classification.tsv.gz", sep="\t", compression="gzip")

    return df_class


def write_stats(df_class):

    stats = df_class.groupby(by="hashID").size().to_dict()
    stats["Total"] = len(df_class)

    with open("stats.yml", "wt") as fout:
        fout.write(yaml.dump(stats))


def parse_arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--hto-umi-count-dir",
        action="store",
        dest="path_hto_umi_count_dir",
        help="path to UMI count outputs generated by CITE-Seq-Count",
        required=True,
    )

    parser.add_argument(
        "--min-count",
        action="store",
        dest="min_count_threshold",
        type=int,
        help="total count for CB less than this threshold will be marked as negative (unreliable observations)",
        default=0,
    )

    parser.add_argument(
        "--mode",
        action="store",
        dest="mode",
        type=int,
        help="processing mode (1=default)",
        default=1,
    )

    # parse arguments
    params = parser.parse_args()

    return params


if __name__ == "__main__":

    params = parse_arguments()

    logger.info("Starting...")

    df_class = hto_demux(
        params.path_hto_umi_count_dir, params.mode, params.min_count_threshold
    )

    logger.info("Writing statistics...")

    write_stats(df_class)

    logger.info("DONE.")

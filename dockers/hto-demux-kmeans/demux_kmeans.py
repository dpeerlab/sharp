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


logger = logging.getLogger("demux_kmeans")

logging.basicConfig(
    level=logging.DEBUG,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler("demux_kmeans.log"),
        logging.StreamHandler(sys.stdout)
    ]
)


def hto_demux(path_dense_count_matrix, path_hto_umi_count_dir):

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

    # convert to numeric cell barcode
    dna3bit = DNA3Bit()
    numeric_barcodes = barcodes.apply(lambda cb: dna3bit.encode(cb))

    df_umi = pd.DataFrame(
        matrix.todense(),
        columns=numeric_barcodes,
        index=features
    ).T

    logger.info(
        "Loaded HTO UMI count matrix ({} x {})".format(
            df_umi.shape[0], df_umi.shape[1]
        )
    )

    # drop the column `unmapped`
    df_umi = df_umi.iloc[:, 0:-1]

    logger.info("Computing centered log-ratio (CLR)...")
    # centered log-ratio (CLR) transformation
    #     	            HTO_301-ACCCACCAGTAAGAC	HTO_302-GGTCGAGAGCATTCA	HTO_303-CTTGCCGCATGTCAT	HTO_304-AAAGCATTCTTCACG
    # 227929296066909	2.609550	0.076485	2.049975	0.137688
    # 164640656084404	2.477301	0.054396	0.046804	3.561632
    # 121748877338358	2.501004	0.091309	0.034176	3.327706
    # 134463437596589	3.060824	2.458869	0.053883
    df_clr = df_umi.apply(lambda row: np.log1p(
        (row + 1) / scipy.stats.mstats.gmean(row + 1)), axis=1)

    # change column name to column index so that we can access by e.g. x[1]
    df_tmp = df_umi
    df_tmp.columns = range(0, 4)

    # for each row (barcode), get the index of the one with the largest UMI count
    ss_umi_largest = df_tmp.idxmax(axis=1)

    def kemans_per_row(row):
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
    df_kmeans = df_clr.apply(lambda row: kemans_per_row(row), axis=1)

    df_kmeans_hotencoded = df_kmeans.apply(
        lambda x: "".join(str(y) for y in x)).to_frame()

    # shorten and replace _ with -
    # ['HTO-301', 'HTO-302', 'HTO-303', 'HTO-304']
    hto_names = list(map(lambda name: name.split(
        "-")[0].replace("_", "-"), df_clr.columns))

    def demux_pass2(cb):

        # index of hto having the largest UMIs: 0, 1, 2, or 3
        idmax = ss_umi_largest[cb]

        # which group belongs to? 0 or 1
        group_id = df_kmeans_hotencoded.loc[cb][0][idmax]

        # how many hto belong that group?
        num_htos = df_kmeans_hotencoded.loc[cb][0].count(group_id)

        # if greater than or equal to two HTOs belong to that group, it means doublet
        # return "Doublet" if doublet, return HTO ID if singlet
        # return "Doublet" if num_htos >= 2 else "Singlet"
        return "Doublet" if num_htos >= 2 else hto_names[idmax]

    df_class = pd.DataFrame(
        list(map(lambda cb: (cb, demux_pass2(cb)), df_kmeans_hotencoded.index))
    )
    df_class.columns = ["CB", "hashID"]
    df_class.set_index("CB", inplace=True)
    df_class

    logger.debug(df_class.groupby(by="hashID").size())

    df_class.to_csv(
        "final-classification.tsv.gz",
        sep="\t",
        compression="gzip"
    )

    df_gene = pd.read_csv(
        path_dense_count_matrix,
        index_col=0
    )

    logger.info(
        "Loaded transcript count matrix ({} x {})".format(
            df_gene.shape[0], df_gene.shape[1]
        )
    )

    df_merged = pd.merge(
        df_gene, df_class,
        left_index=True, right_index=True,
        how="inner"
    )

    logger.info(
        "Merged transcript count matrix with hashtag ({} x {})".format(
            df_merged.shape[0], df_merged.shape[1]
        )
    )

    logger.info("Writing the full dense count matrix with hashtag...")

    df_merged.to_csv(
        "final-matrix.tsv.gz",
        sep="\t",
        compression="gzip"
    )

    return df_class


def write_stats(df_class):

    stats = df_class.groupby(by="hashID").size().to_dict()
    stats["Total"] = len(df_class)

    with open("stats.yml", "wt") as fout:
        fout.write(yaml.dump(stats))


def parse_arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--dense-count-matrix",
        action="store",
        dest="path_dense_count_matrix",
        help="path to scRNA-seq dnese cell-by-gene count matrix file (*.csv)",
        required=True
    )

    parser.add_argument(
        "--hto-umi-count-dir",
        action="store",
        dest="path_hto_umi_count_dir",
        help="path to HTO demux unmapped file (*.csv)",
        required=True
    )

    # parse arguments
    params = parser.parse_args()

    return params


if __name__ == "__main__":

    params = parse_arguments()

    logger.info("Starting...")

    df_class = hto_demux(
        params.path_dense_count_matrix,
        params.path_hto_umi_count_dir
    )

    logger.info("Writing statistics...")
    write_stats(df_class)

    logger.info("DONE.")

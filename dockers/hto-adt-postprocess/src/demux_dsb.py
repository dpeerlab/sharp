#!/usr/bin/env python
# coding: utf-8

import sys
import os
import argparse
import pandas as pd
import numpy as np
from sklearn.mixture import GaussianMixture
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score, davies_bouldin_score
import yaml
import logging
import anndata as ad


logger = logging.getLogger("demux_dsb_kmeans")

logging.basicConfig(
    level=logging.DEBUG,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler("demux_dsb_kmeans.log"),
        logging.StreamHandler(sys.stdout),
    ],
)


def cluster_and_evaluate(data, method="kmeans"):
    """
    Perform clustering and evaluate it using multiple metrics for K-means or GMM.

    Parameters:
    data (np.array): The input data used for clustering
    method (str): Clustering method to use. Either 'kmeans' or 'gmm'. Default is 'kmeans'.

    Returns:
    tuple: A tuple containing:
        - np.array: The cluster labels
        - int: The index of the positive cluster
        - dict: A dictionary containing various goodness of fit metrics
    """
    n_clusters = 2
    if method == "kmeans":
        model = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
        labels = model.fit_predict(data)
        positive_cluster = np.argmax(model.cluster_centers_)

        silhouette = silhouette_score(data, labels)
        davies_bouldin = davies_bouldin_score(data, labels)

        metrics = {"silhouette_score": silhouette, "davies_bouldin_index": davies_bouldin}

    elif method == "gmm":
        model = GaussianMixture(n_components=n_clusters, random_state=42)
        model.fit(data)
        labels = model.predict(data)
        positive_cluster = np.argmax(model.means_)

        bic = model.bic(data)
        log_likelihood = model.score(data) * data.shape[0]

        metrics = {"bic": bic, "log_likelihood": log_likelihood}

    else:
        raise ValueError("Method must be either 'kmeans' or 'gmm'")

    return labels, positive_cluster, metrics

def hto_demux_dsb(
    path_dsb_denoised_adata_dir: str,
    method: str = "kmeans",
):
    """
    Classify HTOs as singlets (assign to HTO), doublets, or negatives based on either a 2-component K-means or GMM,
    and categorize cells based on their HTO classifications.

    Parameters:
    - path_dsb_denoised_adata_dir (str): Path to the DSB denoised anndata directory.
    - method (str): Clustering method to use. Must be either 'gmm' or 'kmeans'. Default is 'kmeans'.

    Returns:
    - AnnData: An AnnData object containing the results of the demultiplexing.
    """
    adata_filtered = ad.read_h5ad(path_dsb_denoised_adata_dir)

    # Check if the dsb_normalized is added in adata_filtered layers
    if "dsb_normalized" in adata_filtered.layers:
        df_umi_dsb_values = adata_filtered.layers["dsb_normalized"]
    else:
        df_umi_dsb_values = adata_filtered.X

    barcodes = adata_filtered.obs_names
    features = adata_filtered.var_names

    df_umi_dsb = pd.DataFrame(
        df_umi_dsb_values,
        columns=features,
        index=barcodes,
    )

    classifications = []
    metrics = {}

    logger.info(f"Running clustering using {method}...")
    for hto in df_umi_dsb.columns:
        data = df_umi_dsb[hto].values.reshape(-1, 1)

        # Perform clustering and evaluation in one step
        labels, positive_cluster, hto_metrics = cluster_and_evaluate(data, method=method)
        metrics[hto] = hto_metrics

        # Classify the points based on the cluster labels
        classifications.append(
            [0 if label != positive_cluster else 1 for label in labels]
        )

    result_df = pd.DataFrame(
        classifications, index=df_umi_dsb.columns, columns=df_umi_dsb.index
    ).T

    # Categorize cells based on their HTO classifications
    def categorize_cell(row):
        positive_htos = row[row == 1].index.tolist()
        if len(positive_htos) == 0:
            return "Negative", None
        elif len(positive_htos) == 1:
            return positive_htos[0][:5], None
        else:
            return "Doublet", ", ".join(positive_htos)

    result_df["hashID"], result_df["Doublet_Info"] = zip(
        *result_df.apply(categorize_cell, axis=1)
    )

    new_df = result_df[["hashID", "Doublet_Info"]]

    logger.info("Classification completed.")

    new_df.to_csv("classification.tsv.gz", sep="\t", compression="gzip")

    # create an anndata object where the denoised data is the X matrix, the barcodes and features are the obs and var names, add the hashID and Doublet_Info as an obs column, and metrics as an uns
    adata_result = ad.AnnData(
        X=df_umi_dsb.values,
        obs=pd.DataFrame(index=df_umi_dsb.index),
        var=pd.DataFrame(index=df_umi_dsb.columns),
    )
    adata_result.obs["hashID"] = result_df["hashID"]
    adata_result.obs["Doublet_Info"] = result_df["Doublet_Info"]
    adata_result.uns["metrics"] = metrics

    return adata_result


def numpy_to_python(obj):
    if isinstance(obj, np.generic):
        return obj.item()
    elif isinstance(obj, dict):
        return {k: numpy_to_python(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [numpy_to_python(i) for i in obj]
    else:
        return obj


def write_stats(result_df, metrics, output_file="stats.yml"):
    stats = result_df.groupby(by="hashID").size().to_dict()
    stats["Total"] = len(result_df)

    # Convert NumPy values to native Python types
    metrics = numpy_to_python(metrics)

    output_dict = {"stats": stats, "metrics": metrics}

    # Write stats and metrics to the YAML file
    with open(output_file, "wt") as fout:
        yaml.dump(output_dict, fout, sort_keys=False, default_flow_style=False)


def parse_arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--dsb-denoised-adata-dir",
        action="store",
        dest="path_dsb_denoised_adata_dir",
        help="path to DSB denoised anndata directory",
        required=True,
    )
    parser.add_argument(
        "--method",
        action="store",
        dest="method",
        help="method used to cluster when demultiplexing",
        default="kmeans",
    )

    parser.add_argument(
        "--output-dir",
        action="store",
        dest="output_dir",
        help="directory to save output files",
        required=True,
    )

    # parse arguments
    params = parser.parse_args()

    return params


if __name__ == "__main__":

    params = parse_arguments()

    logger.info("Starting...")

    adata_result = hto_demux_dsb(
        params.path_dsb_denoised_adata_dir,
        method=params.method,
    )

    os.makedirs(params.output_dir, exist_ok=True)

    stats_file = os.path.join(params.output_dir, "stats.yml")
    write_stats(adata_result.obs, adata_result.uns["metrics"], output_file=stats_file)

    logger.info("Saving AnnData result...")
    adata_output_file = os.path.join(params.output_dir, "demux_result.h5ad")
    adata_result.write(adata_output_file)

    logger.info(f"Results saved to {params.output_dir}")
    logger.info("DONE.")

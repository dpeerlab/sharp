import sys
import argparse
import logging
import anndata as ad

from dsb_algorithm import dsb_adapted

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


def dsb(
    path_adata_filtered_in: str,
    path_adata_raw_in: str,
    path_adata_out: str,
):
    logger.info(f"Loading AnnData {path_adata_filtered_in}...")
    adata_filtered = ad.read_h5ad(path_adata_filtered_in)

    logger.info(f"Loading AnnData {path_adata_raw_in}...")
    adata_raw = ad.read_h5ad(path_adata_raw_in)

    logger.info("Running DSB...")
    dsb_adapted(adata_filtered, adata_raw)

    logger.info(f"Saving AnnData {path_adata_out}...")
    adata_filtered.write(path_adata_out)


def parse_arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--adata-filtered-in",
        action="store",
        dest="path_adata_filtered_in",
        help="path to filtered input AnnData (.h5ad)",
        required=True,
    )

    parser.add_argument(
        "--adata-raw-in",
        action="store",
        dest="path_adata_raw_in",
        help="path to raw input AnnData (.h5ad)",
        required=True,
    )

    parser.add_argument(
        "--adata-out",
        action="store",
        dest="path_adata_out",
        help="path to output AnnData (.h5ad)",
        required=True,
    )

    # parse arguments
    params = parser.parse_args()
    return params


if __name__ == "__main__":

    params = parse_arguments()

    logger.info("Starting...")

    dsb(
        path_adata_filtered_in=params.path_adata_filtered_in,
        path_adata_raw_in=params.path_adata_raw_in,
        path_adata_out=params.path_adata_out,
    )

    logger.info("DONE.")

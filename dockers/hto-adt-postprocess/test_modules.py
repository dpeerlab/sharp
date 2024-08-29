import pytest
import os
import pandas as pd
import anndata as ad

import translate_barcodes
import hto_gex_mapper
import to_adata
import subset_adata


@pytest.fixture
def setup_path():

    path_data = "data"
    path_test_data = "tests"
    yield path_data, path_test_data


def test_hto_gex_mapper_create(setup_path):

    path_data, path_test_data = setup_path

    # $ gunzip -c 3M-february-2018.txt.gz | grep "AAATGGATCGTCGTGA"
    # AAATGGAAGGTCGTGA	AAATGGATCGTCGTGA
    # AAATGGATCGTCGTGA	AAATGGAAGGTCGTGA

    # $ gunzip -c 3M-february-2018.txt.gz | grep "AAATGGATCGTCTTTG"
    # AAATGGAAGGTCTTTG	AAATGGATCGTCTTTG
    # AAATGGATCGTCTTTG	AAATGGAAGGTCTTTG

    # $ gunzip -c 3M-february-2018.txt.gz | grep "TTTGTTGAGTTTCTTC"
    # TTTGTTGAGTTTCTTC	TTTGTTGTCTTTCTTC
    # TTTGTTGTCTTTCTTC	TTTGTTGAGTTTCTTC

    mapper = hto_gex_mapper.create(
        path_10x_whitelist=os.path.join(path_data, "3M-february-2018.txt.gz")
    )

    assert mapper["AAATGGAAGGTCGTGA"] == "AAATGGATCGTCGTGA"
    assert mapper["AAATGGATCGTCTTTG"] == "AAATGGAAGGTCTTTG"
    assert mapper["TTTGTTGAGTTTCTTC"] == "TTTGTTGTCTTTCTTC"


def test_hto_gex_translation_1(setup_path):

    path_data, path_test_data = setup_path

    # $ gunzip -c 3M-february-2018.txt.gz | grep "AAATGGATCGTCGTGA"
    # AAATGGAAGGTCGTGA	AAATGGATCGTCGTGA
    # AAATGGATCGTCGTGA	AAATGGAAGGTCGTGA

    # $ gunzip -c 3M-february-2018.txt.gz | grep "AAATGGATCGTCTTTG"
    # AAATGGAAGGTCTTTG	AAATGGATCGTCTTTG
    # AAATGGATCGTCTTTG	AAATGGAAGGTCTTTG

    barcodes = ["AAATGGAAGGTCGTGA", "AAATGGATCGTCTTTG"]
    df = pd.DataFrame(barcodes).set_index(0)

    translated = translate_barcodes.convert(
        df=df,
        path_hto_gex_mapper=os.path.join(path_data, "10x-hto-gex-mapper.pickle"),
    )

    assert translated.iloc[0].name == "AAATGGATCGTCGTGA"
    assert translated.iloc[1].name == "AAATGGAAGGTCTTTG"


def test_hto_gex_translation_2(setup_path):

    path_data, path_test_data = setup_path

    # $ gunzip -c 3M-february-2018.txt.gz | grep "GCGAGAAGTAGACCGA"
    # GCGAGAACAAGACCGA	GCGAGAAGTAGACCGA
    # GCGAGAAGTAGACCGA	GCGAGAACAAGACCGA

    barcodes = pd.read_csv(
        os.path.join(path_test_data, "barcodes.tsv.gz"),
        sep="\t",
        index_col=0,
        header=None,
        compression="gzip",
    )

    translated = translate_barcodes.convert(
        df=barcodes,
        path_hto_gex_mapper=os.path.join(path_data, "10x-hto-gex-mapper.pickle"),
    )

    for i, barcode in enumerate(barcodes.index):
        if barcode == "GCGAGAAGTAGACCGA":
            assert translated.iloc[i].name == "GCGAGAACAAGACCGA"
            return

    raise Exception("The test file doesn't include the barcode you're testing...")

def get_adata(setup_path, use_acgt=False):
    path_data, path_test_data = setup_path
    path_tag_list = os.path.join(path_test_data, "citeseq", "tag-list.csv")
    path_umi_counts = os.path.join(path_test_data, "citeseq", "umi-counts")

    # to adata
    to_adata.to_adata(sample_name="adata", path_tag_list=path_tag_list, path_umi_counts=path_umi_counts)

    adata = ad.read("adata.h5ad")
    if use_acgt:
        adata.obs_names = adata.obs["barcode_sequence"]
        adata.write("adata.h5ad")

    return adata

def test_to_adata(setup_path):

    adata = get_adata(setup_path)

    assert isinstance(adata, ad.AnnData)
    assert adata.shape[0] > 0
    assert adata.shape[1] > 0
    assert os.path.exists("adata.h5ad")

    os.remove("adata.h5ad")

def test_subset_adata(setup_path):
    path_data, path_test_data = setup_path


    adata = get_adata(setup_path, use_acgt=True)

    path_cb_whitelist = os.path.join(path_test_data, "citeseq", "cb-whitelist.csv")
    cb_whitelist = pd.read_csv(path_cb_whitelist, header=None, index_col=0).index.values

    subset_adata.subset_adata(
        path_adata_in="adata.h5ad",
        path_adata_out="adata.h5ad",
        path_cb_whitelist=path_cb_whitelist,
    )

    adata = ad.read("adata.h5ad")

    os.remove("adata.h5ad")

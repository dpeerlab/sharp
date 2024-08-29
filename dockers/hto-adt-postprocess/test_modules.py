import pytest
import os
import pandas as pd
import anndata as ad

import translate_barcodes
import translate_10x_barcodes
import to_adata
import subset_adata


@pytest.fixture
def base_path(request):
    local = request.config.getoption("--local")
    if local:
        base_path = os.getcwd()
    else:
        base_path = "/opt"
    return base_path


@pytest.fixture
def path_test_data():
    yield "tests"


def test_whitelist_path(base_path):
    """Test whitelist path"""
    path_v3 = translate_10x_barcodes.decide_which_whitelist(
        "test-small-v3", base_path=base_path
    )
    path_v4 = translate_10x_barcodes.decide_which_whitelist(
        "test-small-v4", base_path=base_path
    )

    assert os.path.exists(path_v3)
    assert os.path.exists(path_v4)
    assert path_v3 == os.path.join(base_path, "data/test-small-v3.txt")
    assert path_v4 == os.path.join(base_path, "data/test-small-v4.txt")


def test_translate_barcodes_v3():
    """Test translation of barcodes"""

    # $ gunzip -c data/3M-february-2018.txt.gz | grep "AAACCCAAGAAACACT"
    # $ gunzip -c data/3M-february-2018.txt.gz | grep "AAACCCAAGAAACTGC"

    barcodes = ["AAACCCAAGAAACACT", "AAACCCAAGAAACTGC"]

    translated = translate_barcodes.translate_barcodes(
        barcodes,
        chemistry="test-small-v3",
    )

    assert len(translated) == len(barcodes)
    assert translated[0] == "AAACCCATCAAACACT"
    assert translated[1] == "AAACCCATCAAACTGC"


def test_translate_barcodes_v4():
    """Test translation of barcodes"""

    barcodes = ["AAACCAAAGAACCAGG", "AAACCAAAGAAGCATA"]

    translated = translate_barcodes.translate_barcodes(
        barcodes,
        chemistry="test-small-v4",
    )
    assert len(translated) == len(barcodes)
    assert translated[0] == "AATGAGGTCCATGTCC"
    assert translated[1] == "AATGAGGTCCTGGTAG"


def test_hto_gex_translation_v3():
    """Test V3 translation"""

    # $ gunzip -c 3M-february-2018.txt.gz | grep "AAATGGATCGTCGTGA"
    # AAATGGAAGGTCGTGA	AAATGGATCGTCGTGA
    # AAATGGATCGTCGTGA	AAATGGAAGGTCGTGA

    barcodes = ["AAACCCAAGAAACACT", "AAACCCAAGAAACTGC"]
    df = pd.DataFrame(index=barcodes)

    translated = translate_barcodes.convert(
        df=df,
        chemistry="test-small-v3",
    )

    assert translated.iloc[0].name == "AAACCCATCAAACACT"
    assert translated.iloc[1].name == "AAACCCATCAAACTGC"


def test_hto_gex_translation_v4():
    """Test V4 translation"""

    # $ gunzip -c data/3M-3pgex-may-2023.txt.gz | grep "AAACCAAAGAACCAGG"
    # AAACCAAAGAACCAGG	AATGAGGTCCATGTCC
    # AAACCAAAGAACGGAT	AATGAGGTCTCTAGGG

    barcodes = ["AAACCAAAGAACCAGG", "AAACCAAAGAACGGAT"]
    df = pd.DataFrame(index=barcodes)

    translated = translate_barcodes.convert(
        df=df,
        chemistry="test-small-v4",
    )

    assert translated.iloc[0].name == "AATGAGGTCCATGTCC"
    assert translated.iloc[1].name == "AATGAGGTCTCTAGGG"


def test_hto_gex_translation_large(path_test_data):
    """Test full translation"""

    test_bc = "GCGAGAAGTAGACCGA"

    barcodes = pd.read_csv(
        os.path.join(path_test_data, "barcodes.tsv.gz"),
        sep="\t",
        index_col=0,
        header=None,
        compression="gzip",
    )

    translated = translate_barcodes.convert(df=barcodes, chemistry="10x V3.1 Hashtag")

    assert test_bc in barcodes.index, f"Barcode' {test_bc}' not found in test data..."
    assert translated.index[barcodes.index == test_bc][0] == "GCGAGAACAAGACCGA"


def test_hto_gex_translation_duplicates(path_test_data):
    """Test multiple identical barcodes"""

    barcodes = [
        "AAACCAAAGAACCAGG",
        "AAACCAAAGAACCAGG",
        "AAACCAAAGAACCAGG",
        "AAACCAAAGAACCTAT",
        "AAACCAAAGAACCTAT",
    ]

    translated = translate_barcodes.convert(
        df=pd.DataFrame(index=barcodes),
        chemistry="test-small-v4",
    )

    assert len(translated) == len(
        barcodes
    ), f"Expected {len(barcodes)} barcodes, got {len(translated)}"
    assert translated.iloc[0].name == "AATGAGGTCCATGTCC"
    assert translated.iloc[1].name == "AATGAGGTCCATGTCC"


def test_hto_gex_10x_translation(path_test_data, base_path):

    translate_10x_barcodes.translate(
        path_input=os.path.join(path_test_data, "citeseq", "cb-whitelist-gemx.csv"),
        chemistry="test-small-v4",
        separator=",",
        has_header=False,
        base_path=base_path,
        debug=True,
    )
    assert os.path.exists("translated-barcodes.txt")

    translated_barcodes = pd.read_csv(
        "translated-barcodes.txt", header=None, index_col=0
    )
    assert (
        translated_barcodes.shape[0] == 3
    ), f"Expected 3 barcodes, got {translated_barcodes.shape[0]}"

    os.remove("translated-barcodes.txt")


def get_adata(path_test_data, use_acgt=False):

    path_tag_list = os.path.join(path_test_data, "citeseq", "tag-list.csv")
    path_umi_counts = os.path.join(path_test_data, "citeseq", "umi-counts")

    # to adata
    to_adata.to_adata(
        sample_name="adata",
        path_tag_list=path_tag_list,
        path_umi_counts=path_umi_counts,
    )

    adata = ad.read("adata.h5ad")
    if use_acgt:
        adata.obs_names = adata.obs["barcode_sequence"]
        adata.write("adata.h5ad")

    return adata


def test_to_adata(path_test_data):

    adata = get_adata(path_test_data)

    assert isinstance(adata, ad.AnnData)
    assert adata.shape[0] > 0
    assert adata.shape[1] > 0
    assert os.path.exists("adata.h5ad")

    os.remove("adata.h5ad")


def test_subset_adata(path_test_data):

    adata = get_adata(path_test_data, use_acgt=True)

    path_cb_whitelist = os.path.join(path_test_data, "citeseq", "cb-whitelist.csv")

    subset_adata.subset_adata(
        path_adata_in="adata.h5ad",
        path_adata_out="adata.h5ad",
        path_cb_whitelist=path_cb_whitelist,
        convert=False,
    )

    adata = ad.read("adata.h5ad")
    assert adata is not None

    os.remove("adata.h5ad")

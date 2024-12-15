import pytest
import os
import shutil
import pandas as pd
import anndata as ad

import translate_barcodes
import translate_10x_barcodes
import to_adata
import subset_adata
import fix_barcodes
from tests.utils import get_test_data_path, get_opt_data_path

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


@pytest.fixture
def path_counts_out():
    """Ensure path is empty before and after test."""
    path = "counts"
    if os.path.exists(path):
        shutil.rmtree(path)
    yield path
    if  os.path.exists(path):
        shutil.rmtree(path)


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
    assert path_v3 == get_opt_data_path("test-small-v3.txt")
    assert path_v4 == get_opt_data_path("test-small-v4.txt")

def test_translate_barcodes_v3(base_path):
    """Test translation of barcodes"""
    barcodes = ["AAACCCAAGAAACACT", "AAACCCAAGAAACTGC"]

    translated = translate_barcodes.translate_barcodes(
        barcodes,
        chemistry="test-small-v3",
        base_path=base_path,
    )

    assert len(translated) == len(barcodes)
    assert translated[0] == "AAACCCATCAAACACT"
    assert translated[1] == "AAACCCATCAAACTGC"

def test_translate_barcodes_v4(base_path):
    """Test translation of barcodes"""
    barcodes = ["AAACCAAAGAACCAGG", "AAACCAAAGAAGCATA"]

    translated = translate_barcodes.translate_barcodes(
        barcodes,
        chemistry="test-small-v4",
        base_path=base_path,
    )
    assert len(translated) == len(barcodes)
    assert translated[0] == "AATGAGGTCCATGTCC"
    assert translated[1] == "AATGAGGTCCTGGTAG"

def test_hto_gex_translation_v3(base_path):
    """Test V3 translation"""
    barcodes = ["AAACCCAAGAAACACT", "AAACCCAAGAAACTGC"]
    df = pd.DataFrame(index=barcodes)

    translated = translate_barcodes.convert(
        df=df,
        chemistry="test-small-v3",
        base_path=base_path,
    )

    assert translated.iloc[0].name == "AAACCCATCAAACACT"
    assert translated.iloc[1].name == "AAACCCATCAAACTGC"

def test_hto_gex_translation_v4(base_path):
    """Test V4 translation"""
    barcodes = ["AAACCAAAGAACCAGG", "AAACCAAAGAACGGAT"]
    df = pd.DataFrame(index=barcodes)

    translated = translate_barcodes.convert(
        df=df,
        chemistry="test-small-v4",
        base_path=base_path,
    )

    assert translated.iloc[0].name == "AATGAGGTCCATGTCC"
    assert translated.iloc[1].name == "AATGAGGTCTCTAGGG"

def test_hto_gex_translation_large(base_path):
    """Test full translation"""
    test_bc = "GCGAGAAGTAGACCGA"

    barcodes = pd.read_csv(
        get_test_data_path('tests/tests/barcodes.tsv.gz'),
        sep="\t",
        index_col=0,
        header=None,
        compression="gzip",
    )

    translated = translate_barcodes.convert(
        df=barcodes,
        chemistry="10x V3.1 Hashtag",
        base_path=base_path)

    assert test_bc in barcodes.index, f"Barcode' {test_bc}' not found in test data..."
    assert translated.index[barcodes.index == test_bc][0] == "GCGAGAACAAGACCGA"

def test_hto_gex_translation_duplicates(base_path):
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
        base_path=base_path,
    )

    assert len(translated) == len(barcodes), f"Expected {len(barcodes)} barcodes, got {len(translated)}"
    assert translated.iloc[0].name == "AATGAGGTCCATGTCC"
    assert translated.iloc[1].name == "AATGAGGTCCATGTCC"

def test_hto_gex_10x_translation(path_test_data, base_path):
    translate_10x_barcodes.translate(
        path_input=get_test_data_path('tests/citeseq/cb-whitelist-gemx.csv'),
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
    assert translated_barcodes.shape[0] == 3, f"Expected 3 barcodes, got {translated_barcodes.shape[0]}"

    os.remove("translated-barcodes.txt")

def get_adata(path_test_data, use_acgt=False):
    path_tag_list = get_test_data_path('tests/citeseq/tag-list.csv')
    path_umi_counts = get_test_data_path('tests/citeseq/umi-counts')

    # to adata
    to_adata.to_adata(
        sample_name="adata",
        path_tag_list=path_tag_list,
        path_umi_counts=path_umi_counts,
    )

    adata = ad.read_h5ad("adata.h5ad")
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

    path_cb_whitelist = get_test_data_path('tests/citeseq/cb-whitelist.csv')

    # Non 10x whitelist
    adata = get_adata(path_test_data, use_acgt=True)
    subset_adata.subset_adata(
        path_adata_in="adata.h5ad",
        path_adata_out="adata.h5ad",
        path_cb_whitelist=path_cb_whitelist,
        cb_whitelist_method="not-10x",
        convert=False,
    )

    adata = ad.read_h5ad("adata.h5ad")
    whitelist = pd.read_csv(path_cb_whitelist, header=None).iloc[:,0]
    assert adata is not None
    assert whitelist.isin(adata.obs_names).all()
    os.remove("adata.h5ad")

    # 10x whitelist
    adata = get_adata(path_test_data, use_acgt=True)
    subset_adata.subset_adata(
        path_adata_in="adata.h5ad",
        path_adata_out="adata.h5ad",
        path_cb_whitelist=path_cb_whitelist,
        cb_whitelist_method="10x",
        convert=False,
    )

    adata = ad.read_h5ad("adata.h5ad")
    assert adata.obs_names.str.contains("-1").all
    os.remove("adata.h5ad")



def test_symmetry_of_whitelists(path_test_data):
    """
    Test that all whitelists are symmetrical.
    If this assumption is violated in the future, need to update the translation
    code.
    """

    path_data = get_test_data_path("data")
    whitelists = [f for f in os.listdir(path_data) if f.endswith(".txt.gz")]
    for whitelist in whitelists:
        df = pd.read_csv(os.path.join(path_data, whitelist), names=["gex", "hto"], index_col=None, header=None, sep="\t")
        print(df.head())
        assert all(df.hto.isin(df.gex)), f"ERROR '{whitelist}': Not all HTO barcodes are also GEX barcodes."
        df = df.set_index("gex")
        df.loc[:, "gex_translated"] = df.loc[df.hto].index.values
        assert all(df.hto == df.gex_translated), f"ERROR '{whitelist}': Not all translated GEX barcodes are equivalent to HTO"

@pytest.mark.parametrize("path_counts", [
    get_test_data_path('tests/citeseq/umi-counts'),
    get_test_data_path('tests/citeseq/read-counts')
])
def test_fix_barcodes_not_10x(path_counts_out, path_counts):

    # not 10x
    fix_barcodes.fix_barcodes(
        path_counts_in=path_counts,
        path_counts_out=path_counts_out,
        cb_whitelist_method="not-10x",
    )

    # assertions
    barcodes_in = pd.read_csv(f"{path_counts}/barcodes.tsv.gz", header=None).iloc[:,0]
    barcodes_out = pd.read_csv(f"{path_counts_out}/barcodes.tsv.gz", header=None).iloc[:,0]
    assert not barcodes_in.str.contains("-1").any()
    assert not barcodes_out.str.contains("-1").any()
    assert(barcodes_in == barcodes_out).all()


@pytest.mark.parametrize("path_counts", [
    get_test_data_path('tests/citeseq/umi-counts'),
    get_test_data_path('tests/citeseq/read-counts')
])
def test_fix_barcodes_10x(path_counts_out, path_counts):
    # is 10x
    fix_barcodes.fix_barcodes(
        path_counts_in=path_counts,
        path_counts_out=path_counts_out,
        cb_whitelist_method="10x",
    )

    # assertions
    barcodes_in = pd.read_csv(f"{path_counts}/barcodes.tsv.gz", header=None).iloc[:,0]
    barcodes_out = pd.read_csv(f"{path_counts_out}/barcodes.tsv.gz", header=None).iloc[:,0]
    assert not barcodes_in.str.contains("-1").any()
    assert barcodes_out.str.contains("-1").all()
    assert(barcodes_in + "-1" == barcodes_out).all()

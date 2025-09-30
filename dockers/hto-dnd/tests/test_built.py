"""
Test requirement and correct setup
"""

def test_built():
    import hto
    import anndata
    import scanpy

def test_dnd():
    import hto
    mock = hto.data.generate_hto(n_cells=1000, n_htos=3, noise_level=0.1)
    adata_result = hto.demultiplex(
        adata_hto=mock["filtered"],
        adata_hto_raw=mock["raw"],
        background_method="kmeans-fast",
        add_key_normalise="normalised",
        add_key_denoised="denoised",
        background_version="v2",
        pseudocount=20,
        inplace=False,
    )

    # metrics
    ground_truth = adata_result.obs["ground_truth"]
    predicted = adata_result.obs["hash_id"]

    # evaluate
    precision = hto.metrics.singlet_precision(ground_truth, predicted)
    sensitivity = hto.metrics.singlet_sensitivity(ground_truth, predicted)
    assert precision > 0.9, f"Expected > 0.9, but got {precision:.3f}"
    assert sensitivity > 0.5, f"Expected > 0.5, but got {sensitivity:.3f}"

if __name__ == "__main__":
    try:
        test_built()
        test_dnd()
    except Exception as e:
        print(f"Test failed: {e}")
        exit(1)
    print("All tests passed.")
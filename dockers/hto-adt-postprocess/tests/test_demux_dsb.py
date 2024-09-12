# import pytest
# import numpy as np
# import pandas as pd
# import anndata as ad
# from sklearn.mixture import GaussianMixture
# from sklearn.cluster import KMeans
# from src.demux_dsb import evaluate_clustering, hto_demux_dsb, write_stats

# @pytest.fixture
# def mock_adata():
#     # Create a mock AnnData object
#     np.random.seed(42)
#     n_cells = 100
#     n_htos = 3
#     X = np.random.rand(n_cells, n_htos)
#     obs = pd.DataFrame(index=[f'cell_{i}' for i in range(n_cells)])
#     var = pd.DataFrame(index=[f'HTO_{i}' for i in range(n_htos)])
#     adata = ad.AnnData(X=X, obs=obs, var=var)
#     adata.layers['dsb_normalized'] = X
#     return adata

# def test_evaluate_clustering_kmeans():
#     data = np.array([[1, 2], [1, 4], [1, 0], [4, 2], [4, 4], [4, 0]])
#     results = evaluate_clustering(data, method='kmeans')
#     assert 'silhouette_score' in results
#     assert 'davies_bouldin_index' in results

# def test_evaluate_clustering_gmm():
#     data = np.array([[1, 2], [1, 4], [1, 0], [4, 2], [4, 4], [4, 0]])
#     results = evaluate_clustering(data, method='gmm')
#     assert 'bic' in results
#     assert 'log_likelihood' in results

# def test_hto_demux_dsb(mock_adata, tmp_path):
#     # Save mock AnnData to a temporary file
#     temp_file = tmp_path / "mock_adata.h5ad"
#     mock_adata.write(temp_file)
    
#     # Run hto_demux_dsb
#     result = hto_demux_dsb(str(temp_file), method='kmeans')
    
#     assert isinstance(result, ad.AnnData)
#     assert 'hashID' in result.obs.columns
#     assert 'Doublet_Info' in result.obs.columns
#     assert 'metrics' in result.uns

# def test_write_stats(tmp_path):
#     # Create mock data
#     result_df = pd.DataFrame({
#         'hashID': ['HTO_1', 'HTO_2', 'Negative', 'Doublet'],
#         'Doublet_Info': [None, None, None, 'HTO_1, HTO_2']
#     })
#     metrics = {
#         'HTO_1': {'silhouette_score': 0.8, 'davies_bouldin_index': 0.5},
#         'HTO_2': {'silhouette_score': 0.7, 'davies_bouldin_index': 0.6}
#     }
    
#     # Write stats
#     output_file = tmp_path / "stats.yml"
#     write_stats(result_df, metrics, output_file=str(output_file))
    
#     # Check if file was created
#     assert output_file.exists()

import pytest
import numpy as np
import pandas as pd
import anndata as ad
import yaml
from sklearn.mixture import GaussianMixture
from sklearn.cluster import KMeans
from src.demux_dsb import evaluate_clustering, hto_demux_dsb, write_stats

@pytest.fixture
def mock_dsb_denoised_adata():
    np.random.seed(42)
    n_cells = 1000
    n_htos = 3
    
    # Create bimodal data for each HTO
    data = []
    for _ in range(n_htos):
        background = np.random.normal(loc=0, scale=0.5, size=n_cells // 2)
        signal = np.random.normal(loc=3, scale=0.5, size=n_cells // 2)
        data.append(np.concatenate([background, signal]))
    
    X = np.column_stack(data)
    
    obs = pd.DataFrame(index=[f'cell_{i}' for i in range(n_cells)])
    var = pd.DataFrame(index=[f'HTO_{i}' for i in range(n_htos)])
    
    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.layers['dsb_normalized'] = X
    return adata

def test_evaluate_clustering_kmeans():
    np.random.seed(42)
    data = np.random.gamma(2, 2, size=(1000, 1))
    
    results = evaluate_clustering(data, method='kmeans')
    assert 'silhouette_score' in results
    assert 'davies_bouldin_index' in results
    assert -1 <= results['silhouette_score'] <= 1  # Valid range for silhouette score

def test_evaluate_clustering_gmm():
    np.random.seed(42)
    data = np.random.gamma(2, 2, size=(1000, 1))
    
    results = evaluate_clustering(data, method='gmm')
    assert 'bic' in results
    assert 'log_likelihood' in results
    assert results['log_likelihood'] <= 0  # Log-likelihood should be non-positive

def test_hto_demux_dsb(mock_dsb_denoised_adata, tmp_path):
    temp_file = tmp_path / "mock_dsb_denoised_adata.h5ad"
    mock_dsb_denoised_adata.write(temp_file)
    
    for method in ['kmeans', 'gmm']:
        result = hto_demux_dsb(str(temp_file), method=method)
        
        assert isinstance(result, ad.AnnData)
        assert 'hashID' in result.obs.columns
        assert 'Doublet_Info' in result.obs.columns
        assert 'metrics' in result.uns
        
        # Check if classifications exist
        classifications = result.obs['hashID'].value_counts()
        assert len(classifications) > 0
        assert all(classification in ['Negative', 'Doublet'] or classification.startswith('HTO_') 
                   for classification in classifications.index)
        
        # Check if metrics exist for each HTO
        assert all(f'HTO_{i}' in result.uns['metrics'] for i in range(3))

def test_write_stats(mock_dsb_denoised_adata, tmp_path):
    temp_file = tmp_path / "mock_dsb_denoised_adata.h5ad"
    mock_dsb_denoised_adata.write(temp_file)
    
    result = hto_demux_dsb(str(temp_file), method='kmeans')
    
    output_file = tmp_path / "stats.yml"
    write_stats(result.obs, result.uns['metrics'], output_file=str(output_file))
    
    assert output_file.exists()
    
    # Check content of the YAML file
    with open(output_file, 'r') as f:
        stats = yaml.safe_load(f)
    
    assert 'stats' in stats
    assert 'metrics' in stats
    assert 'Total' in stats['stats']
    assert all(f'HTO_{i}' in stats['metrics'] for i in range(3))

def test_consistent_classification(mock_dsb_denoised_adata, tmp_path):
    temp_file = tmp_path / "mock_dsb_denoised_adata.h5ad"
    mock_dsb_denoised_adata.write(temp_file)
    
    result1 = hto_demux_dsb(str(temp_file), method='kmeans')
    result2 = hto_demux_dsb(str(temp_file), method='kmeans')
    
    pd.testing.assert_series_equal(result1.obs['hashID'], result2.obs['hashID'])
# hto-adt-postprocess

## Update version

* Update config.sh
* Rebuild with `source build.sh`
* Test (see below)
* Push with `source push.sh`

## Unit Tests In Docker

```bash
docker run --rm \
    -v $(pwd)/tests:/opt/tests \
    -v $(pwd)/pytest.ini:/opt/pytest.ini \
    -v $(pwd)/test_modules.py:/opt/test_modules.py \
    sailmskcc/hto-adt-postprocess:0.3.6 \
    pytest -v
```

```
============================= test session starts ==============================
platform linux -- Python 3.9.18, pytest-8.0.0, pluggy-1.4.0 -- /usr/local/bin/python
cachedir: .pytest_cache
rootdir: /opt
configfile: pytest.ini
testpaths: /opt/test_modules.py
collecting ... collected 5 items

test_modules.py::test_hto_gex_mapper_create PASSED                       [ 20%]
test_modules.py::test_hto_gex_translation_1 PASSED                       [ 40%]
test_modules.py::test_hto_gex_translation_2 PASSED                       [ 60%]
test_modules.py::test_to_adata PASSED                                    [ 80%]
test_modules.py::test_subset_adata PASSED                                [100%]

=============================== warnings summary ===============================
test_modules.py::test_to_adata
test_modules.py::test_subset_adata
  /usr/local/lib/python3.9/site-packages/anndata/_core/anndata.py:453: PendingDeprecationWarning: The dtype argument will be deprecated in anndata 0.10.0
    warnings.warn(

-- Docs: https://docs.pytest.org/en/stable/how-to/capture-warnings.html
================== 5 passed, 2 warnings in 122.94s (0:02:02) ===================```

## Run

```bash
docker run -it --rm \
    -v $(pwd)/tests/citeseq:/tests \
    sailmskcc/hto-adt-postprocess:0.3.6
```

### combine.py

```bash
python3 combine.py \
    --dense-count-matrix /data/1187_IL10neg_P163_IGO_09902_8_dense.csv \
    --hto-classification /data/final-classification.tsv.gz
```

### to_adata.py

```bash
python3 to_adata.py \
    --sample test \
    --tag-list /tests/tag-list.csv \
    --umi-counts /tests/umi-counts/ \
    --read-counts /tests/read-counts/
```

docker login docker.io --username sailmskcc # pw: scri1175$$$

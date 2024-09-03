docker run --rm -it \
  -v $(pwd):/opt/project \
  sailmskcc/hto-adt-postprocess:0.3.8 \
  pytest -v /opt/project/tests
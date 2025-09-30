source config.sh

docker run --rm -it \
    --platform linux/amd64 \
    -v $(pwd)/tests:/tests \
    --entrypoint /bin/bash \
    ${image_name}:${version} \
    -c "python /tests/test_built.py"

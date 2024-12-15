# docker-cellranger

Dockerized Cell Ranger v9.0.0

- GEX: https://www.10xgenomics.com/support/software/cell-ranger/downloads

Version 9.0.0 supports GEM-X Flex + CITE, Version 8.0 supports Visium HD.

## License

The code is available to everyone under the standard [MIT license](./LICENSE). However, the code internally uses 10x software, so please make sure that you read and agree to [10x End User Software License](https://www.10xgenomics.com/end-user-software-license-agreement).

## Build Container Image

1. Open to `config.sh`
2. Update the `download_url` (due to expiration of the download link)
3. Run `./build.sh`
4. Run `./push.sh`
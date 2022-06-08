## Prerequisites that should be installed in the docker file.

The prerequisites are installed through [micromamba](https://github.com/mamba-org/micromamba-docker).

- Python:
  - Pandas
  - Biopython

  
## Scripts

- `split_gbk.py` - Split multi-loci GBK file to separate one-locus GBK files.

  
## Build Docker

```
docker buildx create --name mybuilder
docker buildx use mybuilder
docker buildx build --platform=linux/amd64,linux/arm64 -t semenleyn/refprep_split_gbk:latest --push .
```

## Usage

```
docker run \
            --rm \
            -v "$(pwd)/data:/data" \
            -w /data semenleyn/split_gbk \
                                         -g [GBK file] \
                                         -o [Output dicrecory]
```
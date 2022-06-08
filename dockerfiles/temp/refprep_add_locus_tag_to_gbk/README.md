## Prerequisites that should be installed in the docker file.

The prerequisites are installed through [micromamba](https://github.com/mamba-org/micromamba-docker) using [Conda file]('../conda_env/env.yaml').
To reuse layer some prerequisites are installed for related scripts.

- Python:
  - Pandas
  - Biopython

  
## Script

- `add_lt_to_gbk.py` - Adds `/locus_tag` and `/gene` fields to the features in a GBK file.

  
## Build Docker

```
docker buildx create --name mybuilder
docker buildx use mybuilder
docker buildx build --platform=linux/amd64,linux/arm64 -t semenleyn/refprep_add_locus_tag_to_gbk:latest --push .
```

## Usage

```
docker run \
            --rm \
            -v "$(pwd)/data:/data" \
            -w /data semenleyn/refprep_add_locus_tag_to_gbk \
                                                            -g [GBK File] \
                                                            -l [features with locus tags and gene names] \
                                                            -o [Output GBK File]
```
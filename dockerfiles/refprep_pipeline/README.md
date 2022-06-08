# RefPrep docker image

Docker image with tools for running preparation of reference genome for variant calling pipeline.

## Content

- BWA
- Picard
- Samtools
- MUMmer
- Perl
  - JSON module
- Python3
  - Pandas
  - Biopython

## Build

Build using Buildx module: [Docker Buildx](https://docs.docker.com/buildx/working-with-buildx/).
The build is using multi-stage build setup: [Use multi-stage builds](https://docs.docker.com/develop/develop-images/multistage-build/).

```
docker buildx create --name mybuilder
docker buildx use mybuilder
docker buildx build --platform=linux/amd64,linux/arm64 -t semenleyn/refprep_base_conda_env:latest --push .
```


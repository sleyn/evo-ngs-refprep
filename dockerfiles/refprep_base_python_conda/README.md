# Base image for scripts

Used to reduce total images size.

## Prerequisites that should be installed in the docker file.

The prerequisites are installed through [micromamba](https://github.com/mamba-org/micromamba-docker) using [Conda file]('../conda_env/env.yaml').

- Python:
  - Pandas
  - Biopython
  
- Perl
  - JSON


  
## Build Docker

```
docker buildx create --name mybuilder
docker buildx use mybuilder
docker buildx build --platform=linux/amd64,linux/arm64 -t semenleyn/refprep_base_conda_env:latest --push .
```

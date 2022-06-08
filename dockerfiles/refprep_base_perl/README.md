# Base image for Perl scripts

Used to run Perl scripts in pipeline

## Prerequisites that should be installed in the docker file.

Modules:
  - JSON
  
## Build Docker

```
docker buildx create --name mybuilder
docker buildx use mybuilder
docker buildx build --platform=linux/amd64,linux/arm64 -t semenleyn/refprep_base_perl:latest --push .
```

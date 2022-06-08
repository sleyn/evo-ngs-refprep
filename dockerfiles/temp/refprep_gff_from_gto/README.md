## Prerequisites that should be installed in the docker file.

The prerequisites are installed through [micromamba](https://github.com/mamba-org/micromamba-docker).

- Python:
  - Pandas
  - Biopython
  
- Perl:
  - JSON (`conda install -c bioconda perl-json`)
  
## Scripts

- `fna_from_gto.pl` - Convert PATRIC GTO to FNA
- `gff_from_gto_sorted.pl` - Convert PATRIC GTO to GFF
- `add_lt_to_gbk.py` - Adds `/locus_tag` field to the genes in the GBK files.
- `split_gbk.py` - Split multi-loci GBK file to separate one-locus GBK files.

  
## Build Docker

```
docker buildx create --name mybuilder
docker buildx use mybuilder
docker buildx build --platform=linux/amd64,linux/arm64 -t semenleyn/refgenome_fmt_transformers:latest --push .
```

## Usage

```
docker run --rm -v "$(pwd)/data:/data" -w /data semenleyn/refgenome_fmt_transformers [script command]
```
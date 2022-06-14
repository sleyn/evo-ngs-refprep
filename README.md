# Reference prearation workflow for analysis NGS data form   bacterial evolution experiments

The Nextflow workflow designed to generate files required for NGS variant calling of population and clonal data from morbidostat experiments.

## Requirements

The workflow organized using Nextflow. All reuired programs are installed in Docker images.

The pipeline require:

- Java to run Nextflow ( [Java Download](https://www.oracle.com/java/technologies/downloads/) )
- Nextflow ( [Nextflow installation](https://www.nextflow.io/docs/latest/getstarted.html) )
- Docker ( [Docker installation](https://docs.docker.com/get-docker/) )

## Result directory structure

The resulted reference directory should have a following structure:

```
[genome_name]
├── CNOGpro
│   └── [contig].gbk
├── Asample
│   └── vcf
│   │   └── [reactor]A.vcf
│   └── ijump
│   │   └── [reactor]A.txt
│   └── CNOGpro
│       └── CNOGpro_[reactor]A_chr[contig].txt
├── [genome_name].dict
├── [genome_name].fna
├── [genome_name].fna.amb
├── [genome_name].fna.ann
├── [genome_name].fna.bwt
├── [genome_name].fna.fai
├── [genome_name].fna.pac
├── [genome_name].fna.sa
├── [genome_name].gbk
├── [genome_name].gff
├── [genome_name].repeats
├── [genome_name].snpeff_config
├── [genome_name].separator
└── ISTable_processing.txt
```

- `[genome_name].fna` - Nucleotide fasta of genome contigs.
- `[genome_name].gbk` - GeneBank  format file.
- `[genome_name].gff` - GFF3 format file.
- `[genome_name].dict` - PICARD dictionary.
- `[genome_name].fna.fai` - Samtools index.
- `[genome_name].fna.[amb|ann|bwt|pac|sa]` - BWA index files.
- `[genome_name].repeats` - Coordinates of repeat regions for variant filtering.
- `ISTable_processing.txt` - Coordinates of IS elements.
- `[genome_name].snpeff_config` - Generated record for snpEff config file.
- `[genome_name].separator` - File containing information on locus tags structure. Used for CNOGpro **NOTE: Subject for deprication.**
- `CNOGpro` directory - The CNOGpro library that is used for CNV calling can't process multiple contigs. In order to use it we separate original GeneBank file to separate GeneBank file with for each contig.
- `Asample` - Folder with results on variant calling for unevolved `A` samples.

## Usage. The modes

The reference could be prepared with three related workflows. The choise of the workflow is dependent on how reference genome could be abtained.

Modes:

- `gto`
- `rasttk`
- `explicit`

### `gto` mode workflow

The easiest mode - preferred **if a good public genome is available**. It downloads and process all required information from the [PATRIC database](https://www.patricbrc.org/) based on the PATRIC Genome ID (PATRIC genome IDs look like `2345.1` where `2345` is genome id and `1` is genome version).

![gto mode ](/img/diagram_gto.jpeg)
Figure 1. gto mode workflow chart.

Usage:
```
# from the Genomes directory. One level above the [genome_name] directory.
nextflow [path/to/main.nf] --genome [genome_name] --patric_id [patric_id] --mode gto
```

`[genome_name]` - name of the genome. Example: `Escherichia_coli_K12`
`[patric_id]` - PATRIC Genome ID. Example: `511145.12`

### `reasttk` mode workflow

The workflow for genomes annotated by [PATRIC Annotation RASTtk workflow](https://www.patricbrc.org/app/Annotation). This is useful if the genome is not public or not available in PATRIC.

`rasttk` mode requires some preprocessing:

1. Annotate your genome with [PATRIC Annotation RASTtk workflow](https://www.patricbrc.org/app/Annotation).
2. The least automated part is creating a features file `[genome_name].feature`. It is required to add locus tags and gene names to the GeneBank and GFF files.
   
   The `[genome_name].feature` is a table with following columns:
   - `accession` - contig name
   - `start` - start coordinate of gene
   - `end` - end coordinate of gene
   - `strand` - strand where gene is located
   - `refseq_locus_tag` - locus tag
   - `gene` - gene name

A good start to create the table could be `[Genome name].txt` table form the RASTtk output.

3. Download *Nucleotide fasta - .fna*, *GeneBank file - .gbk*, *GFF file - .gff*.
4. Put all required files into common `Genomes` folder to make following structure (**NOTE: You need to put file into `rasttk` folder**):
```
Genomes
└── [genome_name]
    └── rasttk
        ├── [genome_name].fna
        ├── [genome_name].gff
        ├── [genome_name].gbk
        └── [genome_name].feature
```

![gto mode ](/img/diagram_rasttk.jpeg)
Figure 2. rasttk mode workflow chart.

Usage:
```
# from the Genomes directory. One level above the [genome_name] directory.
nextflow [path/to/main.nf] --genome [genome_name] --mode rasttk
```

`[genome_name]` - name of the genome. Example: `Escherichia_coli_K12`

### `explicit` mode workflow

If you already have FNA, GeneBank and GFF files properly formatted you can use `explicit` mode.

The `[genome_name]` directory should be:
```
Genomes
└── [genome_name]
    ├── [genome_name].fna
    ├── [genome_name].gff
    └── [genome_name].gbk
```
![gto mode ](/img/diagram_explicit.jpeg)
Figure 2. explicit mode workflow chart.

Usage:
```
# from the Genomes directory. One level above the [genome_name] directory.
nextflow [path/to/main.nf] --genome [genome_name] --mode explicit
```

## Manual postprocessing

Although workflows do most of required work to prepare reference some manual manupulations will be required.

### 1. Add genome to snpEff annotation database

In order to have annotations for called variants the variant callign pipeline require a snpEff genome database.

To create it we should do following steps:
1.  Create directory with a name of the genome database in `snpEff/data` directory.
2. Put in this folder GBK file renamed it `genes.gbk`
3. In `snpEff.config` located in the snpEff program directory put a description of the genome in the `Non-standard Databases` section generated in `[genome_name].snpeff_config`. Example:
```
# Abaumannii ATCC19606
Abaumannii_ATCC19606.genome: Acenitobacter baumannii ATCC19606
    Abaumannii_ATCC19606.chromosomes : NZ_KL810966.1, NZ_KL810967.1
    Abaumannii_ATCC19606.NZ_KL810966.1.codonTable: Bacterial_and_Plant_Plastid
    Abaumannii_ATCC19606.NZ_KL810967.1.codonTable: Bacterial_and_Plant_Plastid
```
4. run from the SnpEff folder
```
java -jar snpEff.jar build -genbank -v [Database name]
```
5. Finally in the parent `Genomes` directory where all reference folders are located in the `snpEff.txt` file it is a correspondence table of genome foldername and snpEff database name (usually the same). You should add new genome to this table. **NOTE: Subject for deprication.**

| Genome_directory                  | SnpEff                             |
| --------------------------------- | ---------------------------------- |
| Acinetobacter_baumannii_ATCC17978 | Abaumannii_ATCC17978_Assembly_RAST |

### 2. Make a locus tag separator file

This file is used for preparing conddensed CNOGpro output.

`[genome_name].separator` file with all locus tag prefixes and separator between adjacent locus tags. Example:

 | LT      | separator |
 | ------- | --------- |
 | AUO97b_ | 1         |

The separator is a diference between two adjasent locus tags. Example: `AUO97b_00465` and `AUO97b_00466` have separator of one. `BSU_00470` and `BSU_00475` have separator 5.

### 3. Unevolved variants (A-samples)

After unevolved samples will be sequenced and variant calling results files are available they should be added to the reference.

Create `Asample` directory with `vcf` subdirectory for variant calling files, `ijump` for iJump report files and `CNOGpro` for CNOGpro CNV report files. Name files `[reactor]A.vcf`, `[reactor]A.txt` and `CNOGpro_[reactor]A_chr[contig].txt` respectively.

## GFF format details

 - The header `##sequence-region	[contig name]	[contig start coordinate]	[contig end coordinate]` is required for each contig in the reference.
  - Info field should have the following fields:
    - **ID** - Gene ID. For genomes annotated by PATRIC it is PATRIC peg ID.
    - **product** - Gene function.
    - **locus_tag** - Locus tag for the gene.
    - **gene** (optional) - Trivial name.

  Example:

  ```
  ##gff-version 3								
  ##sequence-region	NODE_1_length_3909467_cov_533.478_ID_22129	1	3909467					
  NODE_1_length_3909467_cov_533.478_ID_22129	FIG	CDS	34	336	.	-	1	ID=fig|400667.82.peg.1;product=hypothetical protein;locus_tag=AUO97b_00141
  NODE_1_length_3909467_cov_533.478_ID_22129	FIG	CDS	352	1578	.	-	1	ID=fig|400667.82.peg.2;product=phage replication protein Cri;locus_tag=AUO97b_00142;gene=cri
  NODE_1_length_3909467_cov_533.478_ID_22129	FIG	CDS	1724	2098	.	+	2	ID=fig|400667.82.peg.3;product=helix-turn-helix family protein;locus_tag=AUO97b_00143
  ```
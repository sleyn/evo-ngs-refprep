/*
 * Process GTO.1: Download GTO file
 */

 process DOWNLOAD_GTO {
     tag "${patric_id}"
     container "semenleyn/patric_cli_1.039_ubuntu_20.04:latest"

    input:
      val patric_id
    
    output:
      path "${patric_id}.gto", emit: gto_ch
    
    script:
    """
    p3-gto ${patric_id}
    """
 }

 /*
  * Process GTO.2a: Convert GTO file to FNA
  */

 process GTO_TO_FNA {
    tag "${gto.baseName} gto2fna"
    publishDir "${launchDir}/${params.genome}", mode: 'copy'

    input:
      path gto

    output:
      path "${params.genome}.fna", emit: fna_ch

    script:
    """
    fna_from_gto.pl ${gto} ${params.genome}.fna
    """
 }

 /*
  * Process GTO.2b: Convert GTO file to GFF
  */

process GTO_TO_GFF {
    tag "${gto.baseName} gto2fna"
    publishDir "${launchDir}/${params.genome}", mode: 'copy'

    input:
      path gto

    output:
      path "${params.genome}.gff", emit: gff_ch_raw

    script:
    """
    gff_from_gto_sorted.pl ${gto} ${params.genome}.gff
    """
}

/*
 * Process GTO.2c: Convert GTO file to GBK
 */

process GTO_TO_GBK {
    tag "${gto.baseName} gto2gbk"
    publishDir "${launchDir}/${params.genome}", mode: 'copy'
    container 'semenleyn/patric_cli_1.039_ubuntu_20.04:latest'

    input:
      path gto

    output:
      path "${params.genome}.gbk", emit: gbk_ch_raw

    script:
    """
    rast-export-genome genbank < ${gto} > ${params.genome}.gbk
    """
}

/*
 * Process GTO.3: Generate features file
 */

 process GENERATE_FEATURES {
    tag "${params.patric_id} feature file"
    container "semenleyn/patric_cli_1.039_ubuntu_20.04:latest"

    input:
      val patric_id

    output:
      path "${patric_id}.features", emit: features_ch

    script:
    """
    p3-echo -t genome_id ${patric_id} > ${patric_id}.txt
    p3-get-genome-features --input ${patric_id}.txt \
                           --attr accession,start,end,strand,refseq_locus_tag,gene > ${patric_id}.features
    sed -i.bak 's/feature\\.//g' ${patric_id}.features
    """
 }

/*
 * Process NON-EXPLICIT.1: Add locus tags and gene names to GBK file
 */

process ADD_LT_TO_GBK {
    tag "${gbk.baseName} add lt to gbk"
    publishDir "${launchDir}/${params.genome}", mode: 'copy'

    input:
      path gbk
      path features
    
    output:
      path "${params.genome}.gbk", emit: gbk_ch

    script:
    """
    add_lt_to_gbk.py -g ${gbk} -l ${features} -o ${params.genome}.gbk
    """
}

/*
 * Process NON-EXPLICIT.2: Add locus tags and gene names to GFF file
 */

process ADD_LT_TO_GFF {
    tag "${gff.baseName} add lt to gff"
    publishDir "${launchDir}/${params.genome}", mode: 'copy'

    input:
      path gff
      path fna
      path features
    
    output:
      path "${params.genome}.gff", emit: gff_ch
    
    script:
    """
    add_lt_to_gff.py -g ${gff} -f ${features} -n ${fna} -o ${params.genome}.gff
    """
}

/*
 * Process COMMON.1: Create BWA index
 */

 process PREPARE_BWA_INDEX {
    tag "${fna.baseName} bwa index"
    publishDir "${launchDir}/${params.genome}", mode: 'copy'

    input:
      path fna

    output:
      path "${fna}.amb"
      path "${fna}.ann"
      path "${fna}.bwt"
      path "${fna}.pac"
      path "${fna}.sa"

    script:
    """
    bwa index ${fna} 
    """
 }

 /*
  * Process COMMON.2: Create Samtools index
  */

process PREPARE_SAMTOOLS_INDEX {
    tag "${fna.baseName} samtools index"
    publishDir "${launchDir}/${params.genome}", mode: 'copy'

    input:
      path fna

    output:
      path "${fna}.fai"

    script:
    """
    samtools faidx ${fna}
    """
}

 /*
  * Process COMMON.3: Create Picard dictionary
  */

process PREPARE_PICARD_DICT {
    tag "${fna.baseName} picard dictionary"
    publishDir "${launchDir}/${params.genome}", mode: 'copy'

    input:
      path fna

    output:
      path "${fna.baseName}.dict"

    script:
    """
    java -jar /usr/local/bin/picard.jar CreateSequenceDictionary -R ${fna} -O ${fna.baseName}.dict
    """
}

/*
 * Process COMMON.4: Make IS elements coordinates file
 */

process PREPARE_IS_TABLE {
    tag "${fna.baseName} IS table"
    publishDir "${launchDir}/${params.genome}", mode: 'copy'

    input:
      path fna
    
    output:
      path "ISTable_processing.txt"

    script:
    """
    blastn -query ${fna} -db /data/IS -out IS_blast_out -outfmt 6
    isfinder_db_parcer.py -b IS_blast_out -o .
    """
}

/*
 * Process COMMON.5: Make table with coordinates of repeated sequences
 */

process PREPARE_REPEAT_TABLE {
    tag "${fna.baseName} repeats table"
    publishDir "${launchDir}/${params.genome}", mode: 'copy'

    input:
      path fna

    output:
      path "${fna.baseName}.repeats"

    script:
    """
    nucmer -p ${fna.baseName} -nosimplify --maxmatch ${fna} ${fna}
    show-coords -clrT -I 90 -L 65 ${fna.baseName}.delta > ${fna.baseName}.aln
    make_repeats.pl ${fna.baseName}
    """
}

/*
 * Process COMMON.6: Split GBK file for CNOGpro
 */

 process PREPARE_SPLIT_GBK_CNOGPRO {
    tag "${gbk.baseName} split GBK"
    publishDir "${launchDir}/${params.genome}", mode: 'copy'

    input:
      path gbk

    output:
      path "CNOGpro/*.gbk"

    script:
    """
    split_gbk.py -g ${gbk} -o CNOGpro
    """
 }

 /*
 * Process COMMON.7: Prepare a record for snpEff config file
 */

process PREPARE_SNPEFF_CONFIG {
    tag "${params.genome} snpEff config portion"
    publishDir "${launchDir}/${params.genome}", mode: 'copy'

    input:
      path fna

    output:
      path "${params.genome}.snpeff_config"

    script:
    """
    make_snpeff_config.py -f ${fna} -n ${params.genome} -o ${params.genome}.snpeff_config
    """
}

/*
 * Process 1: Create BWA index
 */

 process PREPARE_BWA_INDEX {
    tag "${fna}_bwa"
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
  * Process 2: Create Samtools index
  */

process PREPARE_SAMTOOLS_INDEX {
    tag "${fna}_samtools"
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
  * Process 3: Create Picard dictionary
  */

process PREPARE_PICARD_DICT {
    tag "${fna}_picard"
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
 * Process 4: Make IS elements coordinates file
 */

process PREPARE_IS_TABLE {
    tag "${fna}_is_table"
    publishDir "${launchDir}/${params.genome}", mode: 'copy'

    input:
      path fna
    
    output:
      path "ISTable_processing.txt"

    script:
    """
    blastn -query ${fna} -db /root/data/IS -out IS_blast_out -outfmt 6
    isfinder_db_parcer.py -b IS_blast_out -o .
    """
}

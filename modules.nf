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
      path fna.baseName

    script:
    """
    java -Dpicard.useLegacyParser=false -jar picard.jar CreateSequenceDictionary -R ${fna} -O ${fna.baseName}
    """
}


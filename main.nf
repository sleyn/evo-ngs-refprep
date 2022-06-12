nextflow.enable.dsl = 2

/*
* Define the default parameters for reference preparation
*/

params.genome = "Escherichia_coli_BW25113"
params.fna = "${launchDir}/${params.genome}/${params.genome}.fna"
params.gbk = "${launchDir}/${params.genome}/${params.genome}.gbk"
params.gff = "${launchDir}/${params.genome}/${params.genome}.gff"

log.info """\
R E F E R E N C E   P R E P A R A T I O N

FNA file: $params.fna
GBK file: $params.gbk
GFF file: $params.gff
"""

/*
 * Import modules
 */

include {
    PREPARE_BWA_INDEX;
    PREPARE_SAMTOOLS_INDEX;
    PREPARE_PICARD_DICT;
    PREPARE_IS_TABLE;
    PREPARE_REPEAT_TABLE;
    PREPARE_SPLIT_GBK_CNOGPRO
} from './modules.nf'

/*
 * Pipeline
 */

workflow {
    fna_ch = Channel.fromPath("$params.fna", checkIfExists: true)
    gbk_ch = Channel.fromPath("$params.gbk", checkIfExists: true)
    gff_ch = Channel.fromPath("$params.gff", checkIfExists: true)

    PREPARE_BWA_INDEX(fna_ch)
    PREPARE_SAMTOOLS_INDEX(fna_ch)
    PREPARE_PICARD_DICT(fna_ch)
    PREPARE_IS_TABLE(fna_ch)
    PREPARE_REPEAT_TABLE(fna_ch)
    PREPARE_SPLIT_GBK_CNOGPRO(gbk_ch)
}
nextflow.enable.dsl = 2

/*
* Define the default parameters for reference preparation
*/

params.genome = "Escherichia_coli_BW25113"
params.fna = "$projectDir/$params.genome/$params.genome\.fna"
params.gbk = "$projectDir/$params.genome/$params.genome\.gbk"
params.gff = "$projectDir/$params.genome/$params.genome\.gff"

log.info """\
R E F E R E N C E   P R E P A R A T I O N
version $manifest.version

FNA file: $params.fna
GBK file: $params.gbk
GFF file: $params.gff
"""

/*
 * Import modules
 */

include {
    PREPARE_BWA_INDEX
} from './modules.nf'
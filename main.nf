nextflow.enable.dsl = 2

/*
 * Define the default parameters for reference preparation
 * based on the selected mode:
 * - rasttk - genome was annotated by RASTtk
 * - gto - genome was downloaded from PATRIC in a GTO format
 * - explicit - all manupulations with genome were done manually and FNA, GBK and GFF files are ready to use
 */

params.genome = "Escherichia_coli_BW25113"
params.mode = 'rasttk'

if( params.mode == 'rasttk' ) {
    /*
     * Reference prepared from a genome annotated by RASTtk
     */

    params.fna = "${launchDir}/${params.genome}/rasttk/${params.genome}.fna"
    params.gbk = "${launchDir}/${params.genome}/rasttk/${params.genome}.gbk"
    params.gff = "${launchDir}/${params.genome}/rasttk/${params.genome}.gff"
    params.feature = "${launchDir}/${params.genome}/rasttk/${params.genome}.feature"

    log.info """\
    R E F E R E N C E   P R E P A R A T I O N

    Reference preparation mode: $params.mode

    FNA file: $params.fna
    GBK file: $params.gbk
    GFF file: $params.gff
    FEATURE file: $params.feature
    """
}else if( params.mode == 'gto' ){
    /*
     * Reference prepared from a genome downloaded from the PATRIC database in GTO format
     */

    params.patric_id = "511145.12" // Using E. coli K-12 substr. MG1655 as default
    // params.feature = "${launchDir}/${params.genome}/gto/${params.genome}.feature"

    log.info """\
    R E F E R E N C E   P R E P A R A T I O N

    Reference preparation mode: $params.mode

    PATRIC ID of the genome: $params.patric_id
    Genome name: $params.genome
    """
}else if( params.mode == 'explicit' ){
    /*
     * Reference prepared from a genome with properly formatted FNA, GFF and GBK files
     */

    params.fna = "${launchDir}/${params.genome}/${params.genome}.fna"
    params.gbk = "${launchDir}/${params.genome}/${params.genome}.gbk"
    params.gff = "${launchDir}/${params.genome}/${params.genome}.gff"

    log.info """\
    R E F E R E N C E   P R E P A R A T I O N

    Reference preparation mode: $params.mode

    FNA file: $params.fna
    GBK file: $params.gbk
    GFF file: $params.gff
    """
}

/*
 * Import modules for GTO mode
 */

if( params.mode == 'gto' ){
    include {
        DOWNLOAD_GTO;
        GTO_TO_FNA;
        GTO_TO_GFF;
        GTO_TO_GBK;
        GENERATE_FEATURES
    } from './modules.nf'
}

if( params.mode != 'explicit' ){
    include {
        ADD_LT_TO_GBK;
        ADD_LT_TO_GFF;
    } from './modules.nf'
}

/*
 * Import common modules
 */

include {
    PREPARE_BWA_INDEX;
    PREPARE_SAMTOOLS_INDEX;
    PREPARE_PICARD_DICT;
    PREPARE_IS_TABLE;
    PREPARE_REPEAT_TABLE;
    PREPARE_SPLIT_GBK_CNOGPRO;
    PREPARE_SNPEFF_CONFIG
} from './modules.nf'

/*
 * Pipeline
 */

workflow {
    if( params.mode == 'gto' ){
        genome_id_ch = Channel.value("${params.patric_id}")
        
        DOWNLOAD_GTO(genome_id_ch)
        gto_ch = DOWNLOAD_GTO.out.gto_ch

        GTO_TO_FNA(gto_ch)
        fna_ch = GTO_TO_FNA.out.fna_ch

        GTO_TO_GFF(gto_ch)
        gff_ch_raw = GTO_TO_GFF.out.gff_ch_raw

        GTO_TO_GBK(gto_ch)
        gbk_ch_raw = GTO_TO_GBK.out.gbk_ch_raw

        GENERATE_FEATURES(genome_id_ch)
        features_ch = GENERATE_FEATURES.out.features_ch
    }else if ( params.mode == 'rasttk' ){
        // Copy FNA file to the main directory
        fna_file = file("${params.fna}")
        fna_file.copyTo("${launchDir}/${params.genome}/${params.genome}.fna")
        fna_ch = Channel.fromPath("${launchDir}/${params.genome}/${params.genome}.fna", checkIfExists: true)

        // Copy GBK file to the main directory
        gbk_file = file("${params.gbk}")
        gbk_file.copyTo("${launchDir}/${params.genome}/${params.genome}.gbk")
        gbk_ch_raw = Channel.fromPath("${launchDir}/${params.genome}/${params.genome}.gbk", checkIfExists: true)

        // Copy GFF file to the main directory
        gff_file = file("${params.gff}")
        gff_file.copyTo("${launchDir}/${params.genome}/${params.genome}.gff")
        gff_ch_raw = Channel.fromPath("${launchDir}/${params.genome}/${params.genome}.gff", checkIfExists: true)

        features_ch = Channel.fromPath("${params.feature}", checkIfExists: true)
    }else{
        fna_ch = Channel.fromPath("${params.fna}", checkIfExists: true)
        gbk_ch = Channel.fromPath("${params.gbk}", checkIfExists: true)
//        gff_ch = Channel.fromPath("${params.gff}", checkIfExists: true)
    }

    if(params.mode != 'explicit'){
        ADD_LT_TO_GBK(gbk_ch_raw, features_ch)
        gbk_ch = ADD_LT_TO_GBK.out.gbk_ch
        ADD_LT_TO_GFF(gff_ch_raw, fna_ch, features_ch)
//        gff_ch = ADD_LT_TO_GFF.out.gff_ch
    }

    PREPARE_BWA_INDEX(fna_ch)
    PREPARE_SAMTOOLS_INDEX(fna_ch)
    PREPARE_PICARD_DICT(fna_ch)
    PREPARE_IS_TABLE(fna_ch)
    PREPARE_REPEAT_TABLE(fna_ch)
    PREPARE_SPLIT_GBK_CNOGPRO(gbk_ch)
    PREPARE_SNPEFF_CONFIG(fna_ch)
}

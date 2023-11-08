// source: https://github.com/nf-core/chipseq/blob/51eba00b32885c4d0bec60db3cb0a45eb61e34c5/subworkflows/local/input_check.nf
//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check.nf'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv
    contrastsheet // file: /path/to/contrastsheet.csv

    main:
    // Validate the samplesheet
    // Validate the contrastsheet
    SAMPLESHEET_CHECK ( 
        samplesheet,
        contrastsheet
        )

    SAMPLESHEET_CHECK.out.run_samplesheet
        .splitCsv( header:true, sep:',', strip:true )
        .map { create_fastq_channel(it) }
        .set { reads }

    emit:
        ch_reads = reads
        versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

def create_fastq_channel(LinkedHashMap row) {
    def meta = [:]
    meta.groupID         = row.groupID
    // meta.libraryID      = row.libraryID
    // meta.sample       = row.sample
    meta.single_end    = row.single_end
    // meta.treat_or_ctrl  = row.treat_or_ctrl

    // add path(s) of the fastq file(s) to the meta map
    def fastq_meta = []
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }
    if (meta.single_end) {
        fastq_meta = [ meta, [ file(row.fastq_1) ] ]
    } else {
        if (!file(row.fastq_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
        }
        fastq_meta = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
    }
    return fastq_meta
}

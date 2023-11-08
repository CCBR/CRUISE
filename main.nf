log.info """\
CRUISE 🛳️
=============
NF version    : $nextflow.version
runName       : $workflow.runName
username      : $workflow.userName
configs       : $workflow.configFiles
profile       : $workflow.profile
cmd line      : $workflow.commandLine
start time    : $workflow.start
projectDir    : $workflow.projectDir
launchDir     : $workflow.launchDir
workDir       : $workflow.workDir
homeDir       : $workflow.homeDir
samplesheet   : ${params.samplesheet}
contrastsheet : ${params.contrastsheet}
"""
.stripIndent()

// SUBMODULES
include { INPUT_CHECK } from './subworkflows/local/input_check.nf'
include { TRIM_COUNT  } from './subworkflows/local/trim_count.nf'

workflow {
    // check samplesheet and contrastsheet
    INPUT_CHECK(file(params.samplesheet),file(params.contrastsheet))

    INPUT_CHECK.out
        .ch_reads
        .set { raw_reads }

    ch_count = params.count_table ? file(params.count_table, checkIfExists: true) : null
    if (!ch_count) { // trim reads and run mageck count
        TRIM_COUNT(raw_reads)
        ch_count = TRIM_COUNT.out.count
    }
}

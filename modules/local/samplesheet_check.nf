// adapted from https://github.com/nf-core/chipseq/blob/51eba00b32885c4d0bec60db3cb0a45eb61e34c5/modules/local/samplesheet_check.nf
process SAMPLESHEET_CHECK {
    tag "$samplesheet"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    path samplesheet
    path contrastsheet

    output:
    path 'valid.samplesheet.csv'            , emit: run_samplesheet
    path 'valid.contrastsheet.csv'        , emit: run_contrastsheet
    path "versions.yml"                                , emit: versions

    script:
    """
    check_samplesheet.py \\
        $samplesheet \\
        $contrastsheet
        
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}

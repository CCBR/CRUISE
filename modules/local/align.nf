
process MAGECK_COUNT {
    label 'mageck'
    container 'quay.io/biocontainers/mageck:0.5.9.5--py39h1f90b4d_3'

    input:
        file(lib)
        val(ids)
        path(fastqs)

    output:
      path("*.count.txt"), emit: count
      path("*.count_normalized.txt"), emit: count_norm
      path("*.countsummary.txt"), emit: count_sum

    script:
    """
    mageck count \\
      -l ${lib} \\
      --fastq ${fastqs} \\
      --sample-label ${ids.join(',')} \\
      -n ${params.exp_name}
    """

    stub:
    """
    for ext in count count_normalized countsummary; do
        touch ${params.exp_name}.\${ext}.txt
    done
    uname -a > output.txt
    echo ${task.container} >> output.txt
    """

}

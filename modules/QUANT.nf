process FEATURE_COUNTS {
    publishDir "${params.outdir}/samples/${pid}/featurecounts", mode: "copy"

    input:
    tuple val(pid), path(bams)
    path gtf

    output:
    path '*.txt'
    file '*.summary'

    script:
    xargs = ""
    if (params.paired) {
      xargs = xargs.concat("-p ")
    }
    """
    featureCounts $xargs \
    -T $task.cpus \
    -t 'exon' \
    -g 'gene_id' \
    -a ${gtf} \
    -o '${pid}_counts.txt' \
    ${bams}
    """
}

process FEATURE_COUNTS_MATRIX {
    publishDir "${params.outdir}/data", mode: "copy"

    input:
    path(counts)

    output:
    path '*.txt'

    script:
    """
    count_matrix.R --format "featurecounts" $counts
    """
}

process ESET {
    publishDir "${params.outdir}/data", mode: "copy"

    input:
    path(matrix)

    output:
    path '*.rds'

    script:
    """
    expression_set.R $matrix
    """
}

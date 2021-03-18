process TRIM_GALORE {
    tag "$pid"
    publishDir "${params.outdir}/samples/${pid}/trimgalore", mode: "copy"

    input:
    tuple val(pid), path(reads)

    output:
    tuple val(pid), path('*.fq.gz')
    path '*.{txt,html,zip}'

    script:
    def data = params.paired ? "--paired ${reads[0]} ${reads[1]}" : "${reads}"
    """
    trim_galore \
    --gzip \
    --fastqc \
    ${data}
    """
}

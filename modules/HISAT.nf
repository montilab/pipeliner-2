process HISAT_INDEX {
    tag "$fasta"
    publishDir "${params.outdir}/hisat", mode: "copy"

    input:
    path fasta
    path gtf 

    output:
    path 'index/*.ht2'

    script:
    """
    mkdir index;
    hisat2-build -p $task.cpus \
                 -f $fasta \
                 index/part
    """
}

process HISAT_MAPPING {
    tag "$pid"
    publishDir "${params.outdir}/samples/${pid}/hisat", mode: "copy"

    input:
    tuple val(pid), path(reads)

    output:
    tuple val(pid), path('*.bam')
    tuple val(pid), path('*.log')

    script:
    def data = params.paired ? "-1 ${reads[0]} -2 ${reads[1]}" : "-U ${reads}"
    """
    hisat2 -p $task.cpus \
           --summary-file '${pid}.log' \
           --new-summary \
           -x ${params.index} \
           ${data} \
           -S '${pid}.sam';
    samtools view -S -b '${pid}.sam';
    samtools sort '${pid}.sam' -o '${pid}.bam'
    """
}

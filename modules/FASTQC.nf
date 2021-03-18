process FASTQC {
    tag "$pid"
    publishDir "${params.outdir}/samples/${pid}/fastqc", mode: "copy"

    input:
    tuple val(pid), path(reads)

    output:
    tuple val(pid), path("*.zip")
    tuple val(pid), path("*.html")

    script:
    def data = params.paired ? "${reads[0]} ${reads[1]}" : "${reads}"
    """
    fastqc -o . -t $task.cpus -q ${data}
    """
}

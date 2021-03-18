process MULTIQC {
    publishDir "${params.outdir}/reports", mode: "copy"

    input:
    path 'data*/*'

    output:
    path '*.html'

    script:
    """
    multiqc ${params.outdir} --force \
    --cl_config "extra_fn_clean_exts: ['_val_1', '_val_2', '.read_distribution', '.read_duplication']"
    """
}

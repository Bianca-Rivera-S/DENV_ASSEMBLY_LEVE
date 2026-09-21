process MULTIQC {
    label 'process_single'
    conda "${projectDir}/environment_illumina.yml"

    publishDir "${params.outdir}/illumina/reports", mode: 'copy'

    input:
    path ('fastp/*')
    path ('mapping/*')
    path ('coverage/*')

    output:
    path "multiqc_report.html", emit: report

    script:
    """
    multiqc . \\
        --filename multiqc_report.html \\
        --title "DENV/ZIKV Illumina Assembly Pipeline Report"
    """
}

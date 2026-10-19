process MULTIQC {
    label 'process_single'
    conda "${projectDir}/environment_illumina.yml"

    publishDir "${params.outdir}/${params.ont ? 'ont' : 'illumina'}/reports", mode: 'copy'

    input:
    path ('fastp/*')
    path ('mapping/*')
    path ('coverage/*')

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_report_data", emit: data

    script:
    def platform = params.ont ? "ONT (Nanopore)" : "Illumina"
    """
    multiqc . \
        --filename multiqc_report.html \
        --title "DENV/ZIKV ${platform} Assembly Pipeline Report"
    """
}


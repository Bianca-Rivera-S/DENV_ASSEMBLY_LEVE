process FASTP {
    tag "$meta.id"
    label 'process_medium'

    publishDir "${params.outdir}/illumina/trimming", mode: 'copy', pattern: "*.fastp.fastq.gz"
    publishDir "${params.outdir}/illumina/qc",       mode: 'copy', pattern: "*.{json,html}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fastp.fastq.gz"), emit: reads
    tuple val(meta), path("*.json")           , emit: json
    tuple val(meta), path("*.html")           , emit: html

    script:
    def prefix   = "${meta.id}"
    def read_in  = meta.single_end ? "--in1 ${reads[0]}" : "--in1 ${reads[0]} --in2 ${reads[1]}"
    def read_out = meta.single_end ? "--out1 ${prefix}_trimmed.fastp.fastq.gz" : "--out1 ${prefix}_1.fastp.fastq.gz --out2 ${prefix}_2.fastp.fastq.gz"

    """
    fastp \
        ${read_in} \
        ${read_out} \
        --qualified_quality_phred 20 \
        --unqualified_percent_limit 30 \
        --length_required 50 \
        --detect_adapter_for_pe \
        --cut_front \
        --cut_tail \
        --thread ${task.cpus} \
        --json ${prefix}.fastp.json \
        --html ${prefix}.fastp.html
    """
}


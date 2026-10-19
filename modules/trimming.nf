process IVAR_TRIM {
    tag "$meta.id"
    label 'process_medium'

    publishDir "${params.outdir}/${params.ont ? 'ont' : 'illumina'}/trimming", mode: 'copy', pattern: "*.{trimmed.sorted.bam,trimmed.sorted.bam.bai}"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta_bed), path(bed)

    output:
    tuple val(meta), path("*.trimmed.sorted.bam"), path("*.trimmed.sorted.bam.bai"), emit: bam
    tuple val(meta), path("*.ivar.log")                                            , emit: log

    script:
    def prefix = "${meta.id}"
    def q_param = params.ont ? 0 : 20
    def m_param = params.ont ? 1 : 30
    """
    ivar trim \\
        -i ${bam} \\
        -b ${bed} \\
        -p ${prefix}.trimmed \\
        -m ${m_param} \\
        -q ${q_param} \\
        -e > ${prefix}.ivar.log

    samtools sort -@ ${task.cpus} -o ${prefix}.trimmed.sorted.bam ${prefix}.trimmed.bam
    samtools index -@ ${task.cpus} ${prefix}.trimmed.sorted.bam
    """
}


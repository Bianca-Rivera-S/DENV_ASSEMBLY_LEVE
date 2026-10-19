process MINIMAP2_ONT {
    tag "$meta.id"
    label 'process_medium'

    publishDir "${params.outdir}/ont/mapping", mode: 'copy'

    input:
    tuple val(meta), path(reads)
    tuple val(ref_meta), path(reference)

    output:
    tuple val(meta), path("${meta.id}.sorted.bam"), path("${meta.id}.sorted.bam.bai"), emit: bam
    tuple val(meta), path("${meta.id}.flagstat"), emit: flagstat

    script:
    def prefix = "${meta.id}"
    """
    minimap2 -ax map-ont -t ${task.cpus} ${reference} ${reads} | \
    samtools view -F 4 -b - | \
    samtools sort -@ ${task.cpus} -o ${prefix}.sorted.bam -
    
    samtools index ${prefix}.sorted.bam
    samtools flagstat ${prefix}.sorted.bam > ${prefix}.flagstat
    """
}


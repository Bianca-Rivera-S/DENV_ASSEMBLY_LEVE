process BWA_MEM {
    tag "$meta.id"
    label 'process_high'
    conda "${projectDir}/environment_illumina.yml"

    publishDir "${params.outdir}/illumina/mapping", mode: 'copy', pattern: "*.{bam,bai,flagstat}"

    input:
    tuple val(meta), path(reads)
    tuple val(meta_ref), path(fasta)

    output:
    tuple val(meta), path("*.sorted.bam"), path("*.sorted.bam.bai"), emit: bam
    tuple val(meta), path("*.flagstat"),                           emit: flagstat

    script:
    def prefix = "${meta.id}"
    """
    if [ ! -f "${fasta}.bwt" ]; then
        bwa index ${fasta}
    fi

    bwa mem -t ${task.cpus} -R '@RG\\tID:${prefix}\\tSM:${prefix}\\tPL:ILLUMINA' ${fasta} ${reads} | \\
    samtools view -@ ${task.cpus} -bS - | \\
    samtools sort -@ ${task.cpus} -o ${prefix}.sorted.bam -

    samtools index -@ ${task.cpus} ${prefix}.sorted.bam
    samtools flagstat -@ ${task.cpus} ${prefix}.sorted.bam > ${prefix}.flagstat
    """
}

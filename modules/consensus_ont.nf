process IVAR_CONSENSUS_ONT {
    tag "$meta.id"
    publishDir "${params.outdir}/ont/consensus", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai), path(ref_fasta)

    output:
    tuple val(meta), path("*.consensus.fa"), emit: consensus

    script:
    """
    # 1. Consenso Estándar (Estricto - Filtra posiciones con profundidad menor a 10x)
    samtools mpileup -aa -A -d 0 -Q 0 --reference ${ref_fasta} ${bam} | \
    ivar consensus -p ${meta.id}.consensus

    # 2. Consenso de Alta Sensibilidad (Baja Cobertura - Recupera bases desde 1x)
    samtools mpileup -aa -A -d 0 -Q 0 --reference ${ref_fasta} ${bam} | \
    ivar consensus -p ${meta.id}.low_coverage.consensus -q 1 -m 1 -t 0
    """
}


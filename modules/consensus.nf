process IVAR_CONSENSUS {
    tag "$meta.id"
    label 'process_medium'
    conda "${projectDir}/environment_illumina.yml"

    publishDir "${params.outdir}/illumina/consensus", mode: 'copy', pattern: "*.fa"
    publishDir "${params.outdir}/illumina/coverage",  mode: 'copy', pattern: "*.{summary.txt,global.dist.txt}"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.consensus.fa"), emit: consensus
    tuple val(meta), path("*.summary.txt") , emit: coverage_summary

    script:
    def prefix = "${meta.id}"
    """
    # 1. Calculate coverage metrics using mosdepth
    mosdepth --threads ${task.cpus} --fast-mode ${prefix} ${bam}

    # 2. Call consensus sequence with iVar
    samtools mpileup -aa -A -d 0 -Q 20 ${bam} | \\
    ivar consensus \\
        -p ${prefix}.consensus \\
        -q 20 \\
        -t 0.6 \\
        -m 10 \\
        -n N

    # Format header to match standard sample ID
    sed -i "1s/>.*/>${prefix}/" ${prefix}.consensus.fa
    """
}


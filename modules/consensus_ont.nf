process IVAR_CONSENSUS_ONT {
    tag "$meta.id"

    publishDir "${params.outdir}/ont/consensus", mode: 'copy', pattern: "*.consensus.fa"
    publishDir "${params.outdir}/ont/coverage", mode: 'copy', pattern: "*.summary.txt"

    input:
    tuple val(meta), path(bam), path(bai), path(reference)

    output:
    tuple val(meta), path("${meta.id}.consensus.fa"), emit: consensus
    tuple val(meta), path("${meta.id}.mosdepth.summary.txt"), emit: coverage_summary

    script:
    """
    # 1. Indizar la referencia si no existe
    if [ ! -f ${reference}.fai ]; then
        samtools faidx ${reference}
    fi

    # 2. Llamada de SNPs únicamente con alta calidad de base (Q15)
    bcftools mpileup -A -d 1000 -q 15 -Q 15 -f ${reference} ${bam} | \
    bcftools call -m --ploidy 1 -V indels -ov -o raw_calls.vcf

    bgzip -f raw_calls.vcf
    bcftools index raw_calls.vcf.gz

    # 3. Crear máscara BED para regiones con cobertura < 20x (enmascaramiento estricto)
    bedtools genomecov -bga -ibam ${bam} | \
    awk '\$4 < 20 {print \$1"\t"\$2"\t"\$3}' > low_coverage.bed

    # 4. Generar consenso aplicando variantes y enmascarando zonas < 20x con N
    bcftools consensus \
        -f ${reference} \
        -m low_coverage.bed \
        raw_calls.vcf.gz > ${meta.id}.consensus.fa

    # 5. Renombrar encabezado FASTA con el ID de la muestra
    sed -i '1s/.*/>${meta.id}/' ${meta.id}.consensus.fa

    # 6. Métrica de cobertura con mosdepth
    mosdepth -n --fast-mode ${meta.id} ${bam}
    """
}


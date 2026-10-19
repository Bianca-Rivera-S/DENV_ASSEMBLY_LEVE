#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process FLYE_DE_NOVO {
    tag "$meta.id"
    publishDir "${params.outdir}/ont/denovo", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_denovo_assembly.fasta"), emit: assembly

    script:
    """
    flye --nano-raw ${reads} --out-dir flye_out --threads ${task.cpus} --gsize 11k
    cp flye_out/assembly.fasta ./${meta.id}_denovo_assembly.fasta
    """
}

process SPADES_DE_NOVO {
    tag "$meta.id"
    publishDir "${params.outdir}/illumina/denovo", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_denovo_assembly.fasta"), emit: assembly

    script:
    """
    spades.py -1 ${reads[0]} -2 ${reads[1]} -o spades_out --threads ${task.cpus}
    cp spades_out/contigs.fasta ./${meta.id}_denovo_assembly.fasta
    """
}


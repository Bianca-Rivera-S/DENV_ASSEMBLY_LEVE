#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
========================================================================================
    IMPORT MODULES
========================================================================================
*/
include { FASTP }              from './modules/qc.nf'
include { UNZIP_ONT }          from './modules/unzip.nf'
include { BWA_MEM }            from './modules/mapping.nf'
include { MINIMAP2_ONT }       from './modules/mapping_ont.nf'
include { IVAR_TRIM }          from './modules/trimming.nf'
include { IVAR_CONSENSUS }     from './modules/consensus.nf'
include { IVAR_CONSENSUS_ONT } from './modules/consensus_ont.nf'
include { MULTIQC }            from './modules/report.nf'

/*
========================================================================================
    PROCESS: AUTOMATIC SEROTYPE DETECTION (ILLUMINA & ONT)
========================================================================================
*/
process DETECT_SEROTYPE {
    tag "$meta.id"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), stdout, emit: serotype

    script:
    """
    if [ ! -f ${projectDir}/references/refs_combined.fasta.bwt ]; then
        cat ${projectDir}/references/*/*.fasta > refs_combined.fasta
        bwa index refs_combined.fasta > /dev/null 2>&1
        REF_DB="refs_combined.fasta"
    else
        REF_DB="${projectDir}/references/refs_combined.fasta"
    fi

    bwa mem -t ${task.cpus} \$REF_DB ${reads[0] ?: reads} 2>/dev/null | \
    samtools view -F 4 -b - | \
    samtools idxstats - | \
    sort -k3 -nr | \
    head -n 1 | \
    cut -f1 > serotype_header.txt

    HEADER=\$(cat serotype_header.txt)

    if [[ "\$HEADER" == *"NC_001477"* ]] || [[ "\$HEADER" == *"DENV1"* ]]; then
        SERO="DENV1"
    elif [[ "\$HEADER" == *"NC_001474"* ]] || [[ "\$HEADER" == *"DENV2"* ]]; then
        SERO="DENV2"
    elif [[ "\$HEADER" == *"NC_001475"* ]] || [[ "\$HEADER" == *"DENV3"* ]]; then
        SERO="DENV3"
    elif [[ "\$HEADER" == *"NC_002640"* ]] || [[ "\$HEADER" == *"DENV4"* ]]; then
        SERO="DENV4"
    else
        SERO="DENV2"
    fi

    echo -n "\$SERO"
    """
}

process DETECT_SEROTYPE_ONT {
    tag "$meta.id"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), stdout, emit: serotype

    script:
    """
    if [ ! -f refs_combined.fasta ]; then
        cat ${projectDir}/references/*/*.fasta > refs_combined.fasta
    fi

    minimap2 -ax map-ont -t ${task.cpus} refs_combined.fasta ${reads} 2>/dev/null | \
    samtools view -F 4 -b - | \
    samtools idxstats - | \
    sort -k3 -nr | \
    head -n 1 | \
    cut -f1 > serotype_header.txt

    HEADER=\$(cat serotype_header.txt)

    if [[ "\$HEADER" == *"NC_001477"* ]] || [[ "\$HEADER" == *"DENV1"* ]]; then
        SERO="DENV1"
    elif [[ "\$HEADER" == *"NC_001474"* ]] || [[ "\$HEADER" == *"DENV2"* ]]; then
        SERO="DENV2"
    elif [[ "\$HEADER" == *"NC_001475"* ]] || [[ "\$HEADER" == *"DENV3"* ]]; then
        SERO="DENV3"
    elif [[ "\$HEADER" == *"NC_002640"* ]] || [[ "\$HEADER" == *"DENV4"* ]]; then
        SERO="DENV4"
    else
        SERO="DENV2"
    fi

    echo -n "\$SERO"
    """
}

/*
========================================================================================
    WORKFLOWS
========================================================================================
*/
workflow ILLUMINA_PIPELINE {
    take:
    ch_samples

    main:
    FASTP(ch_samples)
    DETECT_SEROTYPE(FASTP.out.reads)

    ch_mapped_input = FASTP.out.reads
        .join(DETECT_SEROTYPE.out.serotype)
        .map { meta, reads, detected_serotype ->
            def serotype_clean = detected_serotype.toString().trim()
            def updated_meta   = meta + [ serotype: serotype_clean ]
            def fasta_path     = file("${params.ref_dir}/${serotype_clean}/${serotype_clean}.fasta")
            def ref_meta       = [ id: serotype_clean ]
            return [ updated_meta, reads, ref_meta, fasta_path ]
        }

    BWA_MEM(
        ch_mapped_input.map { meta, reads, ref_meta, fasta -> [ meta, reads ] },
        ch_mapped_input.map { meta, reads, ref_meta, fasta -> [ ref_meta, fasta ] }
    )

    ch_trim_input = BWA_MEM.out.bam.map { meta, bam, bai ->
        def bed_path = file("${params.ref_dir}/${meta.serotype}/${meta.serotype}.bed")
        def bed_meta = [ id: meta.serotype ]
        return [ meta, bam, bai, bed_meta, bed_path ]
    }

    IVAR_TRIM(
        ch_trim_input.map { meta, bam, bai, bed_meta, bed -> [ meta, bam, bai ] },
        ch_trim_input.map { meta, bam, bai, bed_meta, bed -> [ bed_meta, bed ] }
    )

    IVAR_CONSENSUS(IVAR_TRIM.out.bam)

    ch_fastp_reports = FASTP.out.json.map { meta, json -> json }.collect()
    ch_flagstat_logs = BWA_MEM.out.flagstat.map { meta, flag -> flag }.collect()
    ch_coverage_logs = IVAR_CONSENSUS.out.coverage_summary.map { meta, sum -> sum }.collect()

    MULTIQC(ch_fastp_reports, ch_flagstat_logs, ch_coverage_logs)
}

workflow ONT_PIPELINE {
    take:
    ch_compressed

    main:
    UNZIP_ONT(ch_compressed)
    DETECT_SEROTYPE_ONT(UNZIP_ONT.out.fastq)

    ch_mapped_input = UNZIP_ONT.out.fastq
        .join(DETECT_SEROTYPE_ONT.out.serotype)
        .map { meta, reads, detected_serotype ->
            def serotype_clean = detected_serotype.toString().trim()
            def updated_meta   = meta + [ serotype: serotype_clean ]
            def fasta_path     = file("${params.ref_dir}/${serotype_clean}/${serotype_clean}.fasta")
            def ref_meta       = [ id: serotype_clean ]
            return [ updated_meta, reads, ref_meta, fasta_path ]
        }

    // 1. Mapeo ONT con Minimap2
    MINIMAP2_ONT(
        ch_mapped_input.map { meta, reads, ref_meta, fasta -> [ meta, reads ] },
        ch_mapped_input.map { meta, reads, ref_meta, fasta -> [ ref_meta, fasta ] }
    )

    // 2. Vinculación del archivo BED según el serotipo detectado
    ch_trim_input = MINIMAP2_ONT.out.bam.map { meta, bam, bai ->
        def bed_path = file("${params.ref_dir}/${meta.serotype}/${meta.serotype}.bed")
        def bed_meta = [ id: meta.serotype ]
        return [ meta, bam, bai, bed_meta, bed_path ]
    }

    // 3. Recorte de Primers
    IVAR_TRIM(
        ch_trim_input.map { meta, bam, bai, bed_meta, bed -> [ meta, bam, bai ] },
        ch_trim_input.map { meta, bam, bai, bed_meta, bed -> [ bed_meta, bed ] }
    )

    // 4. Preparar entrada para Consenso asociando BAM + referencia FASTA
    ch_consensus_input = IVAR_TRIM.out.bam.map { meta, bam, bai ->
        def fasta_path = file("${params.ref_dir}/${meta.serotype}/${meta.serotype}.fasta")
        return [ meta, bam, bai, fasta_path ]
    }

    // 5. Generación de Consenso
    IVAR_CONSENSUS_ONT(ch_consensus_input)
}

/*
========================================================================================
    MAIN ENTRY POINT
========================================================================================
*/
workflow {
    // Canal Illumina: busca pares R1/R2 o _1/_2
    ch_illumina = Channel.fromFilePairs([
        "${projectDir}/input/illumina/*_{1,2,R1,R2}*.{fastq,fastq.gz,fq,fq.gz}",
        "${projectDir}/input/illumina/**/*_{1,2,R1,R2}*.{fastq,fastq.gz,fq,fq.gz}"
    ])
    .map { pair_id, reads ->
        def clean_id = pair_id.replaceAll(/_S\d+.*$/, '')
        return [ [ id: clean_id, single_end: false ], reads ]
    }

    ILLUMINA_PIPELINE(ch_illumina)

    // Canal Nanopore: soporta .zip, .gz, .fastq, .fq
    ch_ont = Channel.fromPath("${projectDir}/input/ont/*.{zip,gz,fastq,fq}")
        .map { file ->
            def sample_id = file.name.replaceAll(/[-_].*$/, '')
            return [ [ id: sample_id, single_end: true ], file ]
        }

    ONT_PIPELINE(ch_ont)
}


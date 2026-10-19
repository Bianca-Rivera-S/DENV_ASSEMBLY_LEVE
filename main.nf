#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
========================================================================================
    PARÁMETROS POR DEFECTO
========================================================================================
*/
params.illumina = false
params.ont      = false
params.reads    = null
params.serotype = null // Acepta cualquier serotipo/virus (DENV1-4, ZIKA, MAYV, etc.)
params.scheme   = null // Dinámico: No asume esquema de primers por defecto
params.ref_dir  = "references"
params.outdir   = "results"

/*
========================================================================================
    IMPORT MODULES
========================================================================================
*/
include { FASTP }                from './modules/qc.nf'
include { UNZIP_ONT }            from './modules/unzip.nf'
include { BWA_MEM }              from './modules/mapping.nf'
include { MINIMAP2_ONT }         from './modules/mapping_ont.nf'
include { IVAR_TRIM }            from './modules/trimming.nf'
include { IVAR_CONSENSUS }       from './modules/consensus.nf'
include { IVAR_CONSENSUS_ONT }   from './modules/consensus_ont.nf'
include { MULTIQC }              from './modules/report.nf'
include { FLYE_DE_NOVO; SPADES_DE_NOVO } from './modules/denovo.nf'

/*
========================================================================================
    FUNCIÓN AUXILIAR DE RUTAS DE REFERENCIA (DINÁMICA Y COMPATIBLE CON DSL2)
========================================================================================
*/
def getRefFile(serotype, extension) {
    def baseDir = file(params.ref_dir)
    def seroDir = baseDir.resolve(serotype)

    if (!seroDir.exists()) {
        error "Error: No existe el directorio de referencia para '${serotype}' en: ${seroDir.toAbsolutePath()}"
    }

    if (extension == ".bed") {
        // 1. Si se especificó esquema explícito (--scheme 1000 / --scheme 400)
        if (params.scheme) {
            def bedScheme = seroDir.resolve(params.scheme).resolve(serotype + extension)
            if (bedScheme.exists()) return bedScheme
        }

        // 2. Si existe un BED directamente en la raíz de la carpeta del serotipo
        def bedRoot = seroDir.resolve(serotype + extension)
        if (bedRoot.exists()) return bedRoot

        // 3. Auto-detección: Buscar en las subcarpetas de esquemas disponibles usando .find
        def schemeDirs = seroDir.listFiles()?.findAll { it.isDirectory() } ?: []
        def matchBed = schemeDirs.find { dir -> dir.resolve(serotype + extension).exists() }
        if (matchBed) return matchBed.resolve(serotype + extension)

        error "No se encontró archivo BED para '${serotype}' en ${seroDir.toAbsolutePath()}"
    } else {
        // Buscar FASTA en la raíz del directorio del serotipo
        def refFile = seroDir.resolve(serotype + extension)
        if (refFile.exists()) return refFile

        // Buscar FASTA en subcarpetas de esquemas si no está en la raíz
        def schemeDirs = seroDir.listFiles()?.findAll { it.isDirectory() } ?: []
        def matchFasta = schemeDirs.find { dir -> dir.resolve(serotype + extension).exists() }
        if (matchFasta) return matchFasta.resolve(serotype + extension)

        error "No se encontró archivo FASTA para '${serotype}' en ${seroDir.toAbsolutePath()}"
    }
}

/*
========================================================================================
    PROCESSES: AUTOMATIC SEROTYPE DETECTION
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
        cat ${projectDir}/references/*/*.fasta ${projectDir}/references/*/*/*.fasta 2>/dev/null | uniq > refs_combined.fasta
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

    HEADER=\$(cat serotype_header.txt | tr '[:lower:]' '[:upper:]')

    if [[ "\$HEADER" == *"NC_001477"* ]] || [[ "\$HEADER" == *"DENV1"* ]] || [[ "\$HEADER" == *"DENV-1"* ]]; then
        SERO="DENV1"
    elif [[ "\$HEADER" == *"NC_001474"* ]] || [[ "\$HEADER" == *"DENV2"* ]] || [[ "\$HEADER" == *"DENV-2"* ]]; then
        SERO="DENV2"
    elif [[ "\$HEADER" == *"NC_001475"* ]] || [[ "\$HEADER" == *"DENV3"* ]] || [[ "\$HEADER" == *"DENV-3"* ]]; then
        SERO="DENV3"
    elif [[ "\$HEADER" == *"NC_002640"* ]] || [[ "\$HEADER" == *"DENV4"* ]] || [[ "\$HEADER" == *"DENV-4"* ]]; then
        SERO="DENV4"
    elif [[ "\$HEADER" == *"ZIKA"* ]] || [[ "\$HEADER" == *"NC_012532"* ]]; then
        SERO="ZIKA"
    else
        echo "Advertencia: No se identificó serotipo conocido (\$HEADER). Redirigiendo a Ensamblaje De Novo." >&2
        SERO="DE_NOVO"
    fi

    printf "%s" "\$SERO"
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
        cat ${projectDir}/references/*/*.fasta ${projectDir}/references/*/*/*.fasta 2>/dev/null | uniq > refs_combined.fasta
    fi

    minimap2 -ax map-ont -t ${task.cpus} refs_combined.fasta ${reads} 2>/dev/null | \
    samtools view -F 4 -b - | \
    samtools idxstats - | \
    sort -k3 -nr | \
    head -n 1 | \
    cut -f1 > serotype_header.txt

    HEADER=\$(cat serotype_header.txt | tr '[:lower:]' '[:upper:]')

    if [[ "\$HEADER" == *"NC_001477"* ]] || [[ "\$HEADER" == *"DENV1"* ]] || [[ "\$HEADER" == *"DENV-1"* ]]; then
        SERO="DENV1"
    elif [[ "\$HEADER" == *"NC_001474"* ]] || [[ "\$HEADER" == *"DENV2"* ]] || [[ "\$HEADER" == *"DENV-2"* ]]; then
        SERO="DENV2"
    elif [[ "\$HEADER" == *"NC_001475"* ]] || [[ "\$HEADER" == *"DENV3"* ]] || [[ "\$HEADER" == *"DENV-3"* ]]; then
        SERO="DENV3"
    elif [[ "\$HEADER" == *"NC_002640"* ]] || [[ "\$HEADER" == *"DENV4"* ]] || [[ "\$HEADER" == *"DENV-4"* ]]; then
        SERO="DENV4"
    elif [[ "\$HEADER" == *"ZIKA"* ]] || [[ "\$HEADER" == *"NC_012532"* ]]; then
        SERO="ZIKA"
    else
        echo "Advertencia: No se identificó serotipo conocido (\$HEADER). Redirigiendo a Ensamblaje De Novo." >&2
        SERO="DE_NOVO"
    fi

    printf "%s" "\$SERO"
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

    if (params.serotype) {
        ch_serotype = FASTP.out.reads.map { meta, reads -> [ meta, params.serotype ] }
    } else {
        DETECT_SEROTYPE(FASTP.out.reads)
        ch_serotype = DETECT_SEROTYPE.out.serotype
    }

    ch_combined = FASTP.out.reads.join(ch_serotype)

    ch_branched = ch_combined.branch { meta, reads, sero ->
        reference_based: sero.toString().trim() != "DE_NOVO"
        de_novo:         sero.toString().trim() == "DE_NOVO"
    }

    // RUTA A: Mapeo por Referencia
    ch_mapped_input = ch_branched.reference_based.map { meta, reads, detected_serotype ->
        def serotype_clean = detected_serotype.toString().trim()
        def updated_meta   = meta + [ serotype: serotype_clean ]
        def fasta_path     = getRefFile(serotype_clean, ".fasta")
        def ref_meta       = [ id: serotype_clean ]
        return [ updated_meta, reads, ref_meta, fasta_path ]
    }

    BWA_MEM(
        ch_mapped_input.map { meta, reads, ref_meta, fasta -> [ meta, reads ] },
        ch_mapped_input.map { meta, reads, ref_meta, fasta -> [ ref_meta, fasta ] }
    )

    ch_trim_input = BWA_MEM.out.bam.map { meta, bam, bai ->
        def bed_path = getRefFile(meta.serotype, ".bed")
        def bed_meta = [ id: meta.serotype ]
        return [ meta, bam, bai, bed_meta, bed_path ]
    }

    IVAR_TRIM(
        ch_trim_input.map { meta, bam, bai, bed_meta, bed -> [ meta, bam, bai ] },
        ch_trim_input.map { meta, bam, bai, bed_meta, bed -> [ bed_meta, bed ] }
    )

    IVAR_CONSENSUS(IVAR_TRIM.out.bam)

    // RUTA B: Ensamblaje De Novo
    SPADES_DE_NOVO(
        ch_branched.de_novo.map { meta, reads, sero -> [ meta, reads ] }
    )

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

    if (params.serotype) {
        ch_serotype = UNZIP_ONT.out.fastq.map { meta, reads -> [ meta, params.serotype ] }
    } else {
        DETECT_SEROTYPE_ONT(UNZIP_ONT.out.fastq)
        ch_serotype = DETECT_SEROTYPE_ONT.out.serotype
    }

    ch_combined = UNZIP_ONT.out.fastq.join(ch_serotype)

    ch_branched = ch_combined.branch { meta, reads, sero ->
        reference_based: sero.toString().trim() != "DE_NOVO"
        de_novo:         sero.toString().trim() == "DE_NOVO"
    }

    // RUTA A: Mapeo por Referencia
    ch_mapped_input = ch_branched.reference_based.map { meta, reads, detected_serotype ->
        def serotype_clean = detected_serotype.toString().trim()
        def updated_meta   = meta + [ serotype: serotype_clean ]
        def fasta_path     = getRefFile(serotype_clean, ".fasta")
        def ref_meta       = [ id: serotype_clean ]
        return [ updated_meta, reads, ref_meta, fasta_path ]
    }

    MINIMAP2_ONT(
        ch_mapped_input.map { meta, reads, ref_meta, fasta -> [ meta, reads ] },
        ch_mapped_input.map { meta, reads, ref_meta, fasta -> [ ref_meta, fasta ] }
    )

    ch_trim_input = MINIMAP2_ONT.out.bam.map { meta, bam, bai ->
        def bed_path = getRefFile(meta.serotype, ".bed")
        def bed_meta = [ id: meta.serotype ]
        return [ meta, bam, bai, bed_meta, bed_path ]
    }

    IVAR_TRIM(
        ch_trim_input.map { meta, bam, bai, bed_meta, bed -> [ meta, bam, bai ] },
        ch_trim_input.map { meta, bam, bai, bed_meta, bed -> [ bed_meta, bed ] }
    )

    ch_consensus_input = IVAR_TRIM.out.bam.map { meta, bam, bai ->
        def fasta_path = getRefFile(meta.serotype, ".fasta")
        return [ meta, bam, bai, fasta_path ]
    }

    IVAR_CONSENSUS_ONT(ch_consensus_input)

    // RUTA B: Ensamblaje De Novo
    FLYE_DE_NOVO(
        ch_branched.de_novo.map { meta, reads, sero -> [ meta, reads ] }
    )

    ch_fastp_reports = Channel.of(file("NO_FILE_FASTP"))
    ch_flagstat_logs = MINIMAP2_ONT.out.flagstat.map { meta, flag -> flag }.collect()
    ch_coverage_logs = IVAR_CONSENSUS_ONT.out.coverage_summary.map { meta, sum -> sum }.collect()

    MULTIQC(ch_fastp_reports, ch_flagstat_logs, ch_coverage_logs)
}

/*
========================================================================================
    MAIN ENTRY POINT
========================================================================================
*/
workflow {
    if (params.illumina) {
        def read_path = params.reads ?: "${projectDir}/input/illumina/*_R{1,2}_*.fastq.gz"
        ch_illumina = Channel.fromFilePairs(read_path)
            .map { pair_id, reads ->
                def clean_id = pair_id.replaceAll(/_S\d+.*$/, '')
                return [ [ id: clean_id, single_end: false ], reads ]
            }

        ILLUMINA_PIPELINE(ch_illumina)
    }

    if (params.ont) {
        def read_path = params.reads ?: "${projectDir}/input/ont/**/*.{zip,gz,fastq,fq}"
        ch_ont = Channel.fromPath(read_path)
            .map { file ->
                def sample_id = file.parent.name.contains("barcode") ? file.parent.name : file.name.replaceAll(/[-_].*$/, '')
                return [ sample_id, file ]
            }
            .groupTuple()
            .map { sample_id, files ->
                return [ [ id: sample_id, single_end: true ], files ]
            }

        ONT_PIPELINE(ch_ont)
    }

    if (!params.illumina && !params.ont) {
        ch_illumina = Channel.fromFilePairs("${projectDir}/input/illumina/*_R{1,2}_*.fastq.gz")
            .map { pair_id, reads ->
                def clean_id = pair_id.replaceAll(/_S\d+.*$/, '')
                return [ [ id: clean_id, single_end: false ], reads ]
            }

        ILLUMINA_PIPELINE(ch_illumina)

        ch_ont = Channel.fromPath("${projectDir}/input/ont/**/*.{zip,gz,fastq,fq}")
            .map { file ->
                def sample_id = file.parent.name.contains("barcode") ? file.parent.name : file.name.replaceAll(/[-_].*$/, '')
                return [ sample_id, file ]
            }
            .groupTuple()
            .map { sample_id, files ->
                return [ [ id: sample_id, single_end: true ], files ]
            }

        ONT_PIPELINE(ch_ont)
    }
}


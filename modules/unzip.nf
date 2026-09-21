process UNZIP_ONT {
    tag "$meta.id"

    input:
    tuple val(meta), path(infile)

    output:
    tuple val(meta), path("*.fastq"), emit: fastq

    script:
    """
    if [[ "${infile}" == *.zip ]]; then
        # Descomprime el ZIP en un directorio temporal
        unzip -q ${infile} -d unzipped_dir/
        
        # Busca cualquier archivo .fastq o .fq (incluso dentro de subcarpetas)
        FASTQ_FILES=\$(find unzipped_dir/ -type f \\( -name "*.fastq" -o -name "*.fq" -o -name "*.fastq.gz" -o -name "*.fq.gz" \\))
        
        # Concatena todas las lecturas encontradas en un solo FASTQ plano
        zcat \$FASTQ_FILES 2>/dev/null > ${meta.id}.fastq || cat \$FASTQ_FILES > ${meta.id}.fastq
    elif [[ "${infile}" == *.gz ]]; then
        zcat ${infile} > ${meta.id}.fastq
    else
        cp ${infile} ${meta.id}.fastq
    fi
    """
}


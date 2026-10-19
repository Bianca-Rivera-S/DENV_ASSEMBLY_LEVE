process UNZIP_ONT {
    tag "$meta.id"

    input:
    tuple val(meta), path(infiles)

    output:
    tuple val(meta), path("${meta.id}.fastq"), emit: fastq

    script:
    def sample_id = meta.id
    """
    SAMPLE="${sample_id}"
    rm -f \${SAMPLE}.fastq

    for f in ${infiles}; do
        if [[ "\$f" == *.zip ]]; then
            unzip -q "\$f" -d temp_unzip/
            find temp_unzip/ -type f \( -name "*.fastq" -o -name "*.fq" -o -name "*.fastq.gz" -o -name "*.fq.gz" \) | xargs zcat 2>/dev/null >> \${SAMPLE}.fastq || find temp_unzip/ -type f \( -name "*.fastq" -o -name "*.fq" \) | xargs cat >> \${SAMPLE}.fastq
            rm -rf temp_unzip/
        elif [[ "\$f" == *.gz ]]; then
            zcat "\$f" >> \${SAMPLE}.fastq
        else
            cat "\$f" >> \${SAMPLE}.fastq
        fi
    done
    """
}


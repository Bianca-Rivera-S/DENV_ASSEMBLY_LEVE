#!/bin/bash
set -e

# ==========================================
# CONFIGURACIÓN DE MUESTRA Y REFERENCIA
# ==========================================
SAMPLE="barcode01"
REF="references/DENV2/DENV2.fasta"
BED="references/DENV2/DENV2.bed"
THREADS=4

# Directorio de trabajo y salida
OUT="results_manual"
mkdir -p ${OUT}/fastq ${OUT}/mapping ${OUT}/trimming ${OUT}/coverage ${OUT}/consensus

echo "=== 1. CONSOLIDAR ARCHIVOS FASTQ ==="
zcat input/ont/${SAMPLE}/*.fastq.gz > ${OUT}/fastq/${SAMPLE}.fastq
echo "Lecturas descomprimidas exitosamente en: ${OUT}/fastq/${SAMPLE}.fastq"
wc -l ${OUT}/fastq/${SAMPLE}.fastq

echo "=== 2. MAPEO CON MINIMAP2 ==="
minimap2 -ax map-ont -t ${THREADS} ${REF} ${OUT}/fastq/${SAMPLE}.fastq | \
samtools view -F 4 -b - | \
samtools sort -@ ${THREADS} -o ${OUT}/mapping/${SAMPLE}.sorted.bam -

samtools index ${OUT}/mapping/${SAMPLE}.sorted.bam
samtools flagstat ${OUT}/mapping/${SAMPLE}.sorted.bam > ${OUT}/mapping/${SAMPLE}.flagstat
cat ${OUT}/mapping/${SAMPLE}.flagstat | head -n 5

echo "=== 3. RECORTAR PRIMERS CON IVAR TRIM ==="
ivar trim -i ${OUT}/mapping/${SAMPLE}.sorted.bam -b ${BED} -p ${OUT}/trimming/${SAMPLE}.trimmed -m 1 -q 0 -e

samtools sort -@ ${THREADS} -o ${OUT}/trimming/${SAMPLE}.trimmed.sorted.bam ${OUT}/trimming/${SAMPLE}.trimmed.bam
samtools index ${OUT}/trimming/${SAMPLE}.trimmed.sorted.bam

echo "=== 4. COBERTURA CON MOSDEPTH ==="
mosdepth -n --fast-mode ${OUT}/coverage/${SAMPLE} ${OUT}/trimming/${SAMPLE}.trimmed.sorted.bam
cat ${OUT}/coverage/${SAMPLE}.mosdepth.summary.txt

echo "=== 5. GENERACIÓN DE CONSENSO CON IVAR ==="
samtools mpileup -A -d 0 -Q 0 ${OUT}/trimming/${SAMPLE}.trimmed.sorted.bam | \
ivar consensus -p ${OUT}/consensus/${SAMPLE}.consensus -q 10 -m 1 -n N

echo "=== PROCESO COMPLETADO ==="
echo "Primeras líneas del consenso:"
head -n 5 ${OUT}/consensus/${SAMPLE}.consensus.fa


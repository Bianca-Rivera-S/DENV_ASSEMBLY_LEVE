Introduction

DENV_ASSEMBLY_LEVE is a bioinformatics pipeline developed in Nextflow designed to process high-throughput sequencing data (Illumina paired-end and Oxford Nanopore Technologies - ONT) from Dengue virus amplicon schemes.

The workflow accepts raw FASTQ files, performs quality control, adapter/primer trimming, serotype-specific reference mapping, and generates consensus sequences alongside integrated coverage and depth reports.

Tools and Dependencies

This pipeline integrates the following bioinformatic tools:

Nextflow (Workflow orchestration and parallelization)

FastQC / NanoPlot (Read quality control)

Trimmomatic / Chopper (Read filtering and trimming)

BWA-MEM / Minimap2 (Read alignment and reference mapping)

Samtools (BAM processing and depth calculations)

iVar / BCFTools (BED-based primer trimming and variant calling for consensus generation)

MultiQC (Aggregate quality control reporting)


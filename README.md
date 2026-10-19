# DENV_ASSEMBLY_LEVE

## Overview

DENV_ASSEMBLY_LEVE is a bioinformatics pipeline developed using Nextflow for the processing, reference-guided assembly, and genomic characterization of Dengue virus (DENV) sequencing data.

The pipeline is designed to support sequencing data generated using Illumina paired-end and Oxford Nanopore Technologies (ONT), with independent workflows adapted to the characteristics of each sequencing platform.

The main objective is to provide a reproducible, automated, and scalable workflow for generating high-quality DENV consensus genomes from amplicon-based sequencing data.

The pipeline integrates quality control, read filtering, adapter and primer trimming, reference mapping, consensus sequence generation, genome coverage assessment, and quality reporting.

It is designed to support genomic surveillance, molecular epidemiology, phylogenetic analyses, and evolutionary studies of the four Dengue virus serotypes (DENV-1, DENV-2, DENV-3, and DENV-4).

![DENV_ASSEMBLY_LEVE](https://github.com/user-attachments/assets/9963ad92-cdcd-4ed6-9606-78a78549af32)

---

## Main Features

- Processing of Illumina paired-end sequencing data.
- Processing of Oxford Nanopore sequencing data.
- Platform-specific quality control and read filtering.
- Adapter and amplicon primer trimming.
- Reference-guided genome assembly.
- Support for all four DENV serotypes.
- Support for different amplicon sequencing schemes.
- Generation of consensus genome sequences in FASTA format.
- Read alignment and indexed BAM file generation.
- Genome coverage and sequencing depth assessment.
- Identification of genomic regions with insufficient coverage.
- Generation of individual and aggregated quality control reports.
- Reproducible workflow execution using Nextflow.
- Modular architecture for future improvements and additional sequencing technologies.

---

## Pipeline Workflow

DENV_ASSEMBLY_LEVE is organized into independent processing workflows according to the sequencing platform.

### 1. Illumina Workflow

The Illumina workflow processes paired-end FASTQ files through the following steps:

1. Input validation
   - Identification of paired-end FASTQ files.
   - Verification of input file structure.
   - Selection of the corresponding DENV reference genome.

2. Quality control and filtering
   - Assessment of sequencing quality.
   - Removal of sequencing adapters.
   - Filtering of low-quality reads.
   - Generation of quality control reports.

3. Reference mapping
   - Alignment of filtered reads against the selected DENV reference genome.
   - Conversion, sorting, and indexing of alignment files.
   - Generation of BAM and BAI files.

4. Primer trimming
   - Identification of amplicon primer positions using BED files.
   - Removal of primer-derived sequences from mapped reads.
   - Preparation of primer-trimmed alignments for consensus generation.

5. Consensus genome generation
   - Assessment of nucleotide coverage.
   - Generation of a reference-guided consensus sequence.
   - Masking of positions that do not meet the minimum coverage requirements.

6. Genome coverage analysis
   - Calculation of sequencing depth.
   - Assessment of genome breadth of coverage.
   - Identification of low-coverage genomic regions.

7. Quality reporting
   - Generation of sample-level quality metrics.
   - Aggregation of quality control results.
   - Production of reports for downstream interpretation.

### 2. Oxford Nanopore Workflow

The ONT workflow is designed to process long-read FASTQ sequencing data.

The workflow includes:

1. Input validation and reference selection.
2. Read quality filtering.
3. Read length filtering.
4. Reference-guided alignment.
5. Amplicon primer trimming.
6. Consensus sequence generation and polishing.
7. Genome coverage and depth assessment.
8. Quality control reporting.

The ONT workflow is designed to account for the error characteristics of long-read sequencing and to support appropriate consensus generation strategies.

---

## Supported Dengue Virus Serotypes

The pipeline is designed to support the four Dengue virus serotypes:

| Serotype | Identifier |
|----------|------------|
| Dengue virus serotype 1 | DENV1 |
| Dengue virus serotype 2 | DENV2 |
| Dengue virus serotype 3 | DENV3 |
| Dengue virus serotype 4 | DENV4 |

Reference genome sequences are organized separately for each serotype.

The reference sequences and associated primer schemes should be verified before analysis to ensure compatibility with the sequencing protocol.

---

## Bioinformatics Tools

The pipeline is designed to integrate the following tools, depending on the sequencing platform and workflow configuration.

| Tool | Function |
|------|----------|
| Nextflow | Workflow management and reproducibility |
| FastQC | Sequencing quality assessment |
| Fastp | Illumina read filtering and adapter trimming |
| BWA-MEM | Illumina read alignment |
| Minimap2 | Reference alignment, particularly for ONT reads |
| Samtools | BAM processing, sorting, and indexing |
| iVar | Amplicon primer trimming and consensus generation |
| Bcftools | Variant processing and consensus-related operations |
| Bedtools | Genomic interval and coverage operations |
| Mosdepth | Sequencing depth and coverage analysis |
| MultiQC | Aggregated quality control reporting |

Additional tools may be incorporated into the workflow as development progresses.

---

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/Bianca-Rivera-S/DENV_ASSEMBLY_LEVE.git

cd DENV_ASSEMBLY_LEVE
```

### 2. Create the Conda environment

```bash
conda env create -f environment.yml
```

### 3. Activate the environment

```bash
conda activate assembly_env
```

### 4. Verify Nextflow installation

```bash
nextflow -version
```

---

## Input Data

The pipeline accepts sequencing reads in FASTQ format.

### Illumina

Paired-end sequencing files should contain forward and reverse reads.

Example:

```text
input/
└── fastq/
    ├── SAMPLE01_R1.fastq.gz
    ├── SAMPLE01_R2.fastq.gz
    ├── SAMPLE02_R1.fastq.gz
    └── SAMPLE02_R2.fastq.gz
```

### Oxford Nanopore

ONT sequencing data should be provided as FASTQ files.

Example:

```text
input/
└── fastq/
    ├── SAMPLE01.fastq.gz
    └── SAMPLE02.fastq.gz
```

FASTQ files should contain reads that have already undergone basecalling.

---

## Running the Pipeline

The following commands illustrate the intended command-line interface.

Available parameters and execution modes must be verified against the implemented Nextflow workflow.

### 1. Illumina Sequencing

Run the Illumina workflow for DENV-2 using a 1000 bp amplicon scheme:

```bash
nextflow run main.nf \
    --illumina \
    --serotype DENV2 \
    --scheme 1000
```

Run the Illumina workflow using a 400 bp amplicon scheme:

```bash
nextflow run main.nf \
    --illumina \
    --serotype DENV1 \
    --scheme 400
```

### 2. Oxford Nanopore Sequencing

Run the ONT workflow for DENV-3:

```bash
nextflow run main.nf \
    --ont \
    --serotype DENV3 \
    --scheme 1000
```

Run the ONT workflow using a 400 bp amplicon scheme:

```bash
nextflow run main.nf \
    --ont \
    --serotype DENV4 \
    --scheme 400
```

### 3. Automatic Serotype Selection

Automatic serotype identification is a planned functionality intended to allow the pipeline to select the most appropriate DENV reference genome.

Proposed command:

```bash
nextflow run main.nf --illumina
```

This functionality requires an implemented and validated serotype classification strategy.

### 4. Custom Input and Output Directories

Example:

```bash
nextflow run main.nf \
    --illumina \
    --serotype DENV2 \
    --scheme 1000 \
    --reads "input/fastq/*_R{1,2}.fastq.gz" \
    --outdir results
```

### 5. Combined Sequencing Platform Support

The long-term design includes automatic selection of the appropriate workflow based on the sequencing platform.

The objective is to support Illumina and ONT data within a unified Nextflow architecture while maintaining independent processing strategies for each platform.

---

## Output Files

The pipeline is designed to produce the following outputs.

### Consensus Sequences

Consensus genome sequences in FASTA format.

These sequences can be used for downstream analyses, including:

- Multiple sequence alignment.
- Phylogenetic reconstruction.
- Phylogeographic analysis.
- Molecular epidemiology.
- Genomic surveillance.

### Alignment Files

Sorted and indexed BAM files containing sequencing reads mapped against the selected reference genome.

### Primer-Trimmed Alignments

BAM files generated after removal of primer-derived sequences.

### Coverage Reports

Genome coverage reports containing metrics such as:

- Mean sequencing depth.
- Genome breadth of coverage.
- Percentage of reference positions meeting coverage thresholds.
- Regions with insufficient sequencing depth.

### Quality Control Reports

Quality control reports summarizing read filtering, sequencing quality, alignment statistics, and consensus generation metrics.

---

## Expected Output Structure

```text
results/
├── illumina/
│   ├── qc/
│   ├── trimming/
│   ├── mapping/
│   ├── consensus/
│   ├── coverage/
│   └── reports/
│
└── ont/
    ├── qc/
    ├── trimming/
    ├── mapping/
    ├── consensus/
    ├── coverage/
    └── reports/
```

The exact directory structure depends on the implemented Nextflow processes and output configuration.

---

## Quality Considerations

Consensus genome quality depends on several factors, including:

- Initial sequencing quality.
- Sequencing depth.
- Genome coverage uniformity.
- Primer scheme compatibility.
- Reference genome selection.
- Read alignment quality.
- Consensus generation thresholds.

Regions with insufficient sequencing support should be represented using ambiguous nucleotides when appropriate.

A consensus sequence should not be interpreted as a complete high-quality genome solely because a FASTA file was successfully generated.

Genome completeness, sequencing depth, and coverage breadth should be evaluated before downstream analysis.

---

## Reproducibility

DENV_ASSEMBLY_LEVE uses Nextflow to organize bioinformatics processes into modular and reproducible workflows.

Nextflow supports:

- Automated process execution.
- Parallel task processing.
- Dependency management through execution environments.
- Intermediate result caching.
- Workflow resumption.
- Traceability of computational analyses.

Workflow reproducibility also requires documentation of reference genomes, primer schemes, software versions, and analysis parameters.

---

## Development Status

DENV_ASSEMBLY_LEVE is under active development.

The development strategy prioritizes:

1. Implementation and validation of the Illumina workflow.
2. Integration and validation of the Oxford Nanopore workflow.
3. Implementation of sequencing platform detection.
4. Evaluation of automatic DENV serotype identification.
5. Standardization of consensus quality metrics.
6. Improvement of reporting and workflow reproducibility.
7. Validation using sequencing datasets with known characteristics.

Features described as planned should not be considered implemented until they have been integrated and tested.

---

## Applications

The pipeline is intended to facilitate Dengue virus genomic research, including:

- Genomic surveillance of circulating DENV strains.
- Analysis of viral genetic diversity.
- Detection and characterization of circulating viral lineages.
- Reconstruction of evolutionary relationships.
- Investigation of viral introductions and local diversification.
- Generation of genomic datasets for molecular epidemiological studies.

---

## Limitations

DENV_ASSEMBLY_LEVE primarily focuses on reference-guided consensus genome assembly.

The quality of the resulting sequences depends on the sequencing data, reference selection, and compatibility of the amplicon primer scheme.

Automatic serotype detection, de novo assembly, and additional advanced analytical functions require separate implementation and validation.

The pipeline does not replace downstream phylogenetic, phylogeographic, or phylodynamic analyses.

---

## Citation

If you use DENV_ASSEMBLY_LEVE in your research, please cite the repository.

A formal software citation will be provided when a stable version is released.

Repository:

https://github.com/Bianca-Rivera-S/DENV_ASSEMBLY_LEVE

---

## License

License information will be provided in the repository.

---

## Contact

For questions, suggestions, or contributions, please use the GitHub repository issue tracker.

https://github.com/Bianca-Rivera-S/DENV_ASSEMBLY_LEVE


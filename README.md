# GenomeAnonymizer

GenomeAnonymizer is a software for anonymizing Short Read Whole Genome Sequencing (WGS/WES) data by generating a Somatic Tumor Twin (STT) dataset from BAM/CRAM files, by eliminating all germline variation and keeping the somatic variation profile. This ensures privacy preservation while allowing for tumor data analysis allowing researchers and organizations to share and analyze tumor genomic data while protecting individual sensitive information.

## Features

- Anonymizes genomic alignments in BAM/CRAM files
- Output in standard formats (BAM/CRAM) for downstream analysis
- Provides custom classes to use as a library in other Genomic Java applications (e.g., for Pileup analysis of tumor-normal pairs)

## Installation and Usage

GenomeAnonymizer is available as a pre-built container image on DockerHub, making it easy to run with either Docker or Singularity. Alternatively, you can use the JAR file directly if you have Java installed.

### Option 1: Docker Container (Recommended)

#### Prerequisites
- Docker installed on your system

#### Pull and Run
```bash
# Pull the image from DockerHub
docker pull ngaitan55/genomeanonymizer

# Run GenomeAnonymizer
docker run --rm genomeanonymizer:latest java -Xmx16g -jar /GenomeAnonymizer.jar \
    -in $PATH/normal_sample_input.bam \
    -it $PATH/tumor_sample_input.fa \
    -o $PATH/output.bam \
    -r $PATH/reference.fa \
    -t 12
```

### Option 2: Singularity Container

Singularity is commonly used in HPC environments where Docker may not be available.

#### Prerequisites
- Singularity installed on your system

#### Pull and Run
```bash
# Pull the image from DockerHub
singularity pull docker://ngaitan55/genomeanonymizer

# Run GenomeAnonymizer
singularity exec genomeanonymizer_latest.sif java -Xmx16g -jar /GenomeAnonymizer.jar \
    -in $PATH/normal_sample_input.bam \
    -it $PATH/tumor_sample_input.fa \
    -o $PATH/output.bam \
    -r $PATH/reference.fa \
    -t 12
```

### Option 3: Local JAR File

If you prefer to run GenomeAnonymizer directly with Java (useful for development or when containers are not available).

#### Prerequisites
- Java 21 or higher installed on your system
- Gradle for building the project

#### Download and Run
```bash
# Building locally: 
git clone https://github.com/Computational-Genomics-BSC/GenomeAnonymizer.git
cd GenomeAnonymizer
./gradlew build
# The JAR will be available at: build/libs/GenomeAnonymizer-<LATEST_VERSION>.jar
ln -s build/libs/GenomeAnonymizer-<LATEST_VERSION>.jar GenomeAnonymizer.jar

# Run GenomeAnonymizer
java -Xmx16g -jar GenomeAnonymizer.jar \
    -in $PATH/normal_sample_input.bam \
    -it $PATH/tumor_sample_input.fa \
    -o $PATH/output.bam \
    -r $PATH/reference.fa \
    -t 12
```

## Command Line Options

Essential parameters:
- `-in`: Normal sample input BAM/CRAM file
- `-it`: Tumor sample input BAM/CRAM file
- `-o`: Output prefix for the output BAM files
- `-r`: Reference genome FASTA file
- `-t`: Number of threads to use for processing (default: 12)
- `-tmpDir`: Temporary directory for intermediate files, to be used when access to `/tmp` is not available
- `-maxDepth`: Maximum depth threshold for reads in a given genomic region (default: 10000). Use `-1` to disable depth filtering (recommended for gene panels and deep sequencing samples).
- `-h`: Display help information

Advanced parameters:
- `-sampleType`: Type of sample — `WGS` for whole genome sequencing or `gene_panel` for gene panel sequencing (default: `WGS`)
- `-fixVAF`: Correct VAF values at somatic sites after anonymization. Only applies in `gene_panel` mode (flag, default: disabled)
- `-minDepthVAF`: Minimum original tumor depth required to apply VAF correction at a somatic site; sites below this threshold are left uncorrected to avoid false variant calls. Only applies in `gene_panel` mode with `-fixVAF` (default: 20)
- `-minMQ`: Minimum mapping quality for reads to be included in the anonymization process (default: 0; scale: 0–60 PHRED)
- `-includeDuplicates`: Include duplicate reads in the output (flag, default: disabled)
- `-maxReadsInMemory`: Maximum number of reads held in memory at once. Lower this on memory-constrained systems (default: 500000)
- `-s`: Seed for random number generation. Use `-1` for a random seed (default: -1)

### Recommended command for gene panels

When processing gene panel samples, depth filtering should be disabled (`-maxDepth -1`), the sample type set to `gene_panel`, and VAF correction enabled (`-fixVAF`) to preserve accurate allele frequencies after anonymization:

```bash
java -Xmx16g -jar GenomeAnonymizer.jar \
    -in $PATH/normal_sample_input.bam \
    -it $PATH/tumor_sample_input.bam \
    -o $PATH/output \
    -r $PATH/reference.fa \
    -t 12 \
    -sampleType gene_panel \
    -maxDepth -1 \
    -fixVAF
```

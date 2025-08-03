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
# Building locally: ./gradlew build
# The JAR will be available at: build/libs/GenomeAnonymizer-1.0.0.jar

# Run GenomeAnonymizer
java -Xmx16g -jar build/libs/GenomeAnonymizer-1.0.0.jar \
    -in $PATH/normal_sample_input.bam \
    -it $PATH/tumor_sample_input.fa \
    -o $PATH/output.bam \
    -r $PATH/reference.fa \
    -t 12
```

## Command Line Options

Common parameters:
- `-in`: Normal sample input BAM/CRAM file
- `-it`: Tumor sample input BAM/CRAM file
- `-o, --output`: Output prefix for the output BAM files
- `-r, --reference`: Reference genome FASTA file
- `-t, --threads`: Number of threads to use for processing (default is 12)
- `-tmpDir`: Temporary directory for intermediate files, to be used when access to `/tmp` is not available
- `-h`: Display help information
